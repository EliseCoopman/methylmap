from unittest import TestCase, TestSuite
import plotly
import plotly.graph_objects as go
from plotly.colors import DEFAULT_PLOTLY_COLORS as plcolors
import itertools
from argparse import ArgumentParser
import pandas as pd
import os
import gzip
import sys
import numpy as np
import pyranges as pr
from pathlib import Path
import logging
import subprocess

class Region(object):
    def __init__(self, region, expand=False):
        try:
            self.chromosome, interval = region.replace(",", "").split(":")
            try:
                # see if just integer chromosomes are used
                self.chromosome = int(self.chromosome)
            except ValueError:
                pass
            self.begin, self.end = [int(i) for i in interval.split("-")]
        except ValueError:
            sys.exit(
                "\n\nERROR: Window (-w/--window) inproperly formatted, "
                "an example of accepted formats is:\n'chr5:150200605-150423790'\n\n"
            )
        if expand:
            self.begin = self.begin - int(expand)
            self.end = self.end + int(expand)
        self.start = self.begin  # start is now just an alias for begin because I tend to forget
        self.size = self.end - self.begin
        if self.size < 0:
            sys.exit(
                "\n\nERROR: Window (-w/--window) inproperly formatted, "
                "begin of the interval has to be smaller than end\n\n"
            )
        self.string = f"{self.chromosome}_{self.begin}_{self.end}"
        self.fmt = f"{self.chromosome}:{self.begin}-{self.end}"


class Transcript(object):
    def __init__(self, transcript, gene, exon_tuples, strand):
        self.transcript = transcript
        self.gene = gene
        self.exon_tuples = list(exon_tuples)
        self.strand = strand
        self.marker = "triangle-up" if self.strand == "+" else "triangle-down"
        self.begin = min(list(itertools.chain.from_iterable(self.exon_tuples)))
        self.end = max(list(itertools.chain.from_iterable(self.exon_tuples)))
        self.color = ""


def main():
    args = get_args()
    meth_browser(
        files=args.files,
        table=args.table,
        outtable=args.outtable,
        outfig=args.outfig,
        names=args.names,
        window=args.window,
        expand=args.expand,
        gff=args.gff,
        groups=args.groups,
        simplify=args.simplify,
    )


def get_args():
    parser = ArgumentParser(description="Create heatmap of methylation frequencies.")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument(
        "-f",
        "--files",
        nargs="+",
        help="Haplotype specific methylation frequency files: nanopolish calculate_methylation_frequency.py output or BAM/CRAM files",
    )
    action.add_argument("-t", "--table", help="methfrequencytable or overviewtable")
    # chr17:44345246-44353106
    parser.add_argument("-w", "--window", help="Region to visualise.")
    parser.add_argument("-n", "--names", nargs="*", default=[])
    parser.add_argument("--gff", "--gtf", help="gtf/gff3 file")
    parser.add_argument(
        "--expand",
        help="number of base pairs to expand the window with in both directions",
        type=int,
        default=0,
    )
    parser.add_argument("--outtable", help="File to write the frequencies table to.")
    parser.add_argument("--outfig", help="File to write output heatmap to.")
    parser.add_argument("--groups", nargs="*", help="List of category per sample")
    parser.add_argument("-s", "--simplify", action="store_true")  # default: False
    args = parser.parse_args()
    if args.files:
        if len(args.names) == 0:
            for b in args.files:
                c = b.rsplit(".", 2)[0]
                args.names.append(os.path.basename(c))
        else:
            if len(args.files) != len(args.names):
                sys.exit(
                    f"ERROR: expecting the same number of input files [{len(args.files)}] and names [{len(args.names)}]"
                )
    return args


def meth_browser(
    files,
    table,
    outfig="browser.html",
    outtable=False,
    names=False,
    window=False,
    expand=False,
    gff=False,
    groups=False,
    simplify=False,
):
    if window:
        window = Region(window, expand)

    num_col = 2 if gff else 1  # number of subplots (columns) needed
    subplots = create_subplots(num_col)

    # frequencies table with all meth frequencies of all samples
    meth_data, window = read_mods(
        files, table, names, window, groups, gff, outtable
    ) 
    fig = plot_methylation(subplots, meth_data, num_col)

    if gff:
        annotation_traces = gff_annotation(gff, window, simplify)
        for annot_trace in annotation_traces:
            fig.append_trace(trace=annot_trace, row=1, col=1)
        fig.update_xaxes(title_text="", showticklabels=False, zeroline=False, row=1, col=1)
        fig.update_yaxes(title_text="", showticklabels=True, zeroline=False, row=1, col=1)
    with open(outfig, "w") as f:
        f.write(fig.to_html())


def read_mods(files, table, names, window, groups, gff, outtable):
    """
    converts a file from nanopolish to a pandas dataframe
    input can be from calculate_methylation_frequency
    which will return a dataframe with 'chromosome', 'pos', 'methylated_frequency'
    """
    if files:
        file_type = file_sniffer(files[0])
    elif table:
        file_type = file_sniffer(table)
    try:
        if file_type == "nanopolish_calc_meth_freq":
            return parse_nanopolish(files, table, names, window, groups, outtable)
        elif file_type == "methfrequencytable":
            return parse_methfrequencytable(table, names, window, groups, gff, outtable)
        elif file_type == "overviewtable_nanopolishfiles":
            return parse_nanopolish(files, table, names, window, groups, outtable)
        # elif file_type in ["cram", "bam"]:
        #     return parse_cram(filename, file_type, name, window)
    except Exception as e:
        logging.error("Error processing input file(s).")
        logging.error(e, exc_info=True)
        sys.stderr.write("\n\n\nError processing input file(s)!\n")
        raise


def parse_overviewtable(table):
    overviewtable = pd.read_table(table).sort_values(["group", "name"])
    files = overviewtable["path"].tolist()
    names = overviewtable["name"].tolist()
    return files, names


def parse_nanopolish(files, table, names, window, groups, outtable):
    """
    input files = calculate_methylation_frequency.py output files
    """
    if table:
        files, names = parse_overviewtable(table)

    dfs = []
    for file, name in zip(files, names):
        if window:
            if not Path(file + ".tbi").is_file():
                try:
                    with open(file) as f:
                        subprocess.Popen(["tabix", "-S1", "-s1", "-b2", "-e3", file], stdout=f)
                except FileNotFoundError as e:
                    logging.error("Error when making a .tbi file.")
                    logging.error(e, exc_info=True)
                    sys.stderr.write("\n\nERROR when making a .tbi file.\n")
                    sys.stderr.write("Is tabix installed and on the PATH?.")
                    raise
            try:
                logging.info(f"Reading {file} using a tabix stream.")
                tabix_stream = subprocess.Popen(
                    ["tabix", file, window.fmt],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
            except FileNotFoundError as e:
                logging.error("Error when opening a tabix stream.")
                logging.error(e, exc_info=True)
                sys.stderr.write("\n\nERROR when opening a tabix stream.\n")
                sys.stderr.write("Is tabix installed and on the PATH?")
                raise
            header = gzip.open(file, "rt").readline().rstrip().split("\t")
            table = pd.read_csv(tabix_stream.stdout, sep="\t", header=None, names=header)
            logging.info("Read the file in a dataframe.")
        else:
            table = pd.read_csv(file, sep="\t")  # telkens hele file inlezen
            logging.info("Read the file in a dataframe.")
        table.drop(["group_sequence","called_sites_methylated","num_motifs_in_group","called_sites"],axis=1,inplace=True)
        table["position"] = (table["start"] + table["end"]) / 2
        table = table.set_index("position").drop(["chromosome", "start", "end"], axis=1)
        dfs.append(table.rename(columns={"methylated_frequency": name}))
    methfrequencytable = dfs[0].join(dfs[1:], how="outer")
    methfreqtable = methfrequencytable.sort_values("position", ascending=False)
    # output is an meth frequency table with position as index and for each sample a column with all the methylation frequencies

    if files:
        if groups:
            headerlist = list(methfreqtable.columns.values)
            if len(headerlist) == len(groups):
                res = zip(headerlist, groups)
                output = sorted(list(res), key=lambda x: x[1])
                orderedlist = [i[0] for i in output]
                methfreqtable = methfreqtable.reindex(columns=orderedlist)
            else:
                sys.exit(
                    f"ERRORwhen matching --groups with samples, is length of --groups list ({len(groups)}) matching with number of sample files?")
    methfreqtable.to_csv(outtable, sep="\t", na_rep=np.NaN, header=True)
    return methfreqtable, window


def parse_methfrequencytable(table, names, window, groups, gff, outtable):
    table = pd.read_table(table).sort_values("position", ascending=False)
    if (
        window
    ):  ### WDC look at between https://pandas.pydata.org/docs/reference/api/pandas.Series.between.html
        methfreqtable = methfreqtable[
            (methfreqtable["position"] >= window.begin) & (methfreqtable["position"] <= window.end)
        ]
    else:
        if gff:
            ### WDC would suggest to make a chrom column required
            if "chrom" in methfreqtable.columns:
                ### WDC I would check if methfreqtable["chrom"].unique() has length of 1
                if (methfreqtable["chrom"] == methfreqtable["chrom"][0]).all():
                    chrom = methfreqtable.iloc[0, methfreqtable.columns.get_loc("chrom")]
                    if not chrom.startswith("chr"):  ### WDC not all chromosomes have 'chr'
                        logging.error(
                            "Error when extracting window out of methfreqtable positions."
                        )
                        sys.stderr.write(
                            "\n\nError when extracting window out of methfreqtable positions.\n"
                        )
                        sys.stderr.write(
                            "Is position column of methfreqtable in format chr1:123456?"
                        )
                else:
                    ### WDC probably have to stop the script if this happens
                    logging.error("Error when extracting window out of methfreqtable positions.")
                    sys.stderr.write(
                        "\n\nError when extracting window out of methfreqtable positions.\n"
                    )
                    sys.stderr.write("Window over different chromosomes is not possible.")
                begin = float(methfreqtable.iloc[numberofpositions, 0])
                end = float(methfreqtable.iloc[0, 0])
                window = Region(f"{chrom}:{round(begin)}-{round(end)}")
            
    methfreqtable.drop(["chrom"], axis=1, inplace=True).set_index("position", inplace=True)

    if names:
        methfreqtable.columns = names

    if groups:
        try:
            headerlist = list(methfreqtable.columns.values)
            res = zip(headerlist, groups)
            output = sorted(list(res), key=lambda x: x[-1])
            orderedlist = [i[0] for i in output]
            methfreqtable = methfreqtable.reindex(columns=orderedlist)
        ### WDC I don't think you could get a FileNotFoundError here?
        except FileNotFoundError as e:
            logging.error("Error when matching --groups with samples.")
            logging.error(e, exc_info=True)
            sys.stderr.write("\n\nError when matching --groups with samples.\n")
            sys.stderr.write(
                f"Is length of the --groups list ({len(groups)}) matching with number of sample files?"
            )
            raise

    methfreqtable.to_csv(outtable, sep="\t", na_rep=np.NaN, header=True)
    return methfreqtable, window


# def parse_cram(files, filetype, names, window): #####nog niet werkende   #modbam2bed #modbamtools #mbtools
#     import pysam
#     mode = 'rc' if filetype == 'cram' else 'rb'
#     cram = pysam.AlignmentFile(filename, mode)
#     data = []
#     start_stops = []
#     for read in cram.fetch(reference=str(window.chromosome), start=window.begin, end=window.end):
#         if not read.is_supplementary and not read.is_secondary:
#             start_stops.append((read.query_name, read.reference_start, read.reference_end))
#             mod_positions = get_modified_reference_positions(read)
#             if mod_positions:
#                 data.extend(mod_positions)
#     df = pd.DataFrame(data, columns=['read_name', 'strand', 'pos', 'quality', 'mod']) \
#         .astype(dtype={'mod': 'category', 'quality': 'float'}) \
#         .sort_values(['read_name', 'pos'])

#     return (table=sub_df, data_type="ont-cram",
#                          name=f"{name}_{mod}",
#                          called_sites=len(sub_df),
#                          start_end_table=pd.DataFrame(start_stops,
#                          columns=['read_name', 'posmin', 'posmax'])
#                          .set_index('read_name')
#                          )
#             for mod, sub_df in df.groupby('mod')


def gff_annotation(gff, window, simplify):
    result = []
    annotation = annotation_transcripts(gff, window, simplify)
    for x_pos, transcript in enumerate(annotation):
        line = make_per_gene_annot_line_trace(transcript, window, x_pos)
        exons = [
            make_per_exon_arrow_trace(transcript, begin, end, x_pos)
            for begin, end in transcript.exon_tuples
            if window.begin < begin and window.end > end
        ]
        result.extend([line, *exons])
    return result


def annot_file_type(annot_file):
    """
    Figure out type of annotation file.
    """
    if annot_file.endswith((".gtf", ".gtf.gz")):
        return "gtf"
    elif annot_file.endswith((".gff", ".gff.gz", ".gff2", ".gff2.gz", ".gff3", ".gff3.gz")):
        return "gff"
    else:
        sys.exit(
            "ERROR: unrecognized extension of the annotation file.\n"
            "Supported are gtf, gtf.gz, gff, gff.gz, gff2, gff2.gz, gff3 and gff3.gz"
        )


def transcripts_in_window(df, window, feature="transcript"):
    """
    Return the transcript names for which
    either the end or the begin of an exon is within the window
    """
    return df.loc[
        df["begin"].between(window.begin, window.end) | df["end"].between(window.begin, window.end),
        feature,
    ].unique()


def assign_colors_to_genes(transcripts):
    genes = set([t.gene for t in transcripts])
    colordict = {g: c for g, c in zip(genes, plcolors * 100)}
    for t in transcripts:
        t.color = colordict[t.gene]


def parse_annotation(gff, window):
    """
    Parse the gff and select the relevant window as determined by the window
    """
    type = annot_file_type(gff)
    logging.info(f"Parsing {type} file...")
    if Path(gff + ".tbi").is_file():
        logging.info(f"Reading {gff} using a tabix stream.")
        try:
            tabix_stream = subprocess.Popen(
                ["tabix", gff, window.fmt],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        except FileNotFoundError as e:
            logging.error("Error when opening a tabix stream.")
            logging.error(e, exc_info=True)
            sys.stderr.write("\n\nERROR when opening a tabix stream.\n")
            sys.stderr.write("Is tabix installed and on the PATH?.")
            raise
    else:
        try:
            with open(gff) as f:
                ### WDC Does this attempt to overwrite the gff?
                ### WDC is there even a stdout with this command?
                subprocess.Popen(
                    ["tabix", "-s1", "-b4", "-e5", gff], stdout=f
                )  # possible problem: -S, -s, -b, -e zijn anders in andere gff files
                # WDC "tabix ,"-p", gff" -> dit werkt niet
        except FileNotFoundError as e:
            logging.error("Error when making a .tbi file.")
            logging.error(e, exc_info=True)
            sys.stderr.write("\n\nERROR when making a .tbi file.\n")
            sys.stderr.write("Is tabix installed and on the PATH?.")
            raise
        return parse_annotation(gff, window, simplify=False)
    annotationfile = pd.read_csv(
        tabix_stream.stdout,
        sep="\t",
        header=None,
        names=[
            "chromosome",
            "source",
            "feature_type",
            "begin",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ],
    )
    return annotationfile


def annotation_transcripts(gff, window, simplify):
    annotationfile = parse_annotation(gff, window)
    gene = []
    transcript = []
    annotationfile["attributes"] = annotationfile.attributes.str.split(";")

    for x in annotationfile["attributes"]:
        transcript_idfound = False
        for i in x:
            if i.startswith("gene_name"):
                gene.append(i.split("=")[1])
            if i.startswith("transcript_id"):
                transcript_idfound = True
                transcript.append(i.split("=")[1])
        if not transcript_idfound:
            transcript.append("")

    annotationfile["gene"] = gene
    annotationfile["transcript"] = transcript

    annotationfile.drop(
        ["source", "feature_type", "score", "phase", "attributes"], axis=1, inplace=True
    )

    if simplify:
        annotationfile.drop_duplicates(subset=["chromosome", "begin", "end", "gene"], inplace=True)
        res = []
        for g in transcripts_in_window(annotationfile, window, feature="gene"):
            gtable = annotationfile.loc[annotationfile["gene"] == g]
            if len(gtable):
                res.append(
                    Transcript(
                        transcript=gtable["gene"].tolist()[0],
                        gene=gtable["gene"].tolist()[0],
                        exon_tuples=gtable.loc[:, ["begin", "end"]]
                        .sort_values("begin")
                        .itertuples(index=False, name=None),
                        strand=gtable["strand"].tolist()[0],
                    )
                )
        sys.stderr.write(f"Found {len(res)} gene(s) in the window.\n")
        logging.info(f"Found {len(res)} gene(s) in the window.\n")
    else:
        res = []
        for t in transcripts_in_window(annotationfile, window, feature="transcript"):
            tr = annotationfile.loc[annotationfile["transcript"] == t]
            res.append(
                Transcript(
                    transcript=t,
                    gene=tr["gene"].tolist()[0],
                    exon_tuples=tr.loc[:, ["begin", "end"]]
                    .sort_values("begin")
                    .itertuples(index=False, name=None),
                    strand=tr["strand"].tolist()[0],
                )
            )
        sys.stderr.write(f"Found {len(res)} transcript(s) in the window.\n")
        logging.info(f"Found {len(res)} transcript(s) in the window.\n")
    assign_colors_to_genes(res)
    return res


def make_per_gene_annot_line_trace(transcript, window, x_pos):
    """Generate a line trace for the gene
    Trace can get limited by the window sizes
    """
    return go.Scatter(
        y=[max(transcript.begin, window.begin), min(transcript.end, window.end)],
        x=[x_pos, x_pos],
        mode="lines",
        line=dict(width=2, color=transcript.color),
        name=transcript.transcript,
        text=transcript.gene,
        hoverinfo="text",
        showlegend=False,
    )


def make_per_exon_arrow_trace(transcript, begin, end, x_pos):
    """Generate a line+marker trace for the exon
    The shape is an arrow, as defined by the strand in transcript.marker
    """
    return go.Scatter(
        y=[begin, end],
        x=[x_pos, x_pos],
        mode="lines+markers",
        line=dict(width=8, color=transcript.color),
        name=transcript.transcript,
        text=transcript.gene,
        hoverinfo="text",
        showlegend=False,
        marker=dict(symbol=transcript.marker, size=8),
    )


def file_sniffer(filename):
    """
    Takes in a filename and tries to guess the input file type
    """
    if not Path(filename).is_file():
        sys.exit(f"\n\nERROR: File {filename} does not exist, please check the path!\n")
    if is_bam_file(filename):  # input BAM
        return "bam"
    if is_cram_file(filename):  # input CRAM
        return "cram"
    # input: calculate_methylation_frequency.py output (.tsv or .tsv.gz) OR own methfrequencytable (.tsv or .tsv.gz)
    if is_gz_file(filename):
        import gzip

        header = gzip.open(filename, "rt").readline()
    else:
        header = open(filename, "r").readline()

    if "methylated_frequency" in header:
        # calculate_methylation_frequency.py output as input
        return "nanopolish_calc_meth_freq"
    if "path" in header:  # overviewtable
        df = pd.read_table(filename)
        if df["path"].iloc[0].endswith(".tsv"):  ###only checks first file in overviewtable
            return "overviewtable_nanopolishfiles"
        elif df["path"].iloc[0].endswith(".tsv.gz"):  ###only checks first file in overviewtable
            return "overviewtable_nanopolishfiles"
        elif df["path"].iloc[0].endswith(".bam"):  ####only checks first file in overviewtable
            return "bam"
        elif df["path"].iloc[0].endswith(".cram"):  ####only checks first file in overviewtable
            return "cram"
    if header.startswith("chromosome"):
        ####own methfrequencytable as input: needs first column header to be "chromosome" (nog niet de ideale oplossing)
        return "methfrequencytable"
    sys.exit(f"\n\n\nInput file {filename} not recognized!\n")


def is_gz_file(filepath):
    import binascii

    with open(filepath, "rb") as test_f:
        return binascii.hexlify(test_f.read(2)) == b"1f8b"


def is_cram_file(filepath):
    with open(filepath, "rb") as test_f:
        return test_f.read(4) == b"CRAM"


def is_bam_file(filepath):
    import gzip

    try:
        with gzip.open(filepath) as test_f:
            return test_f.read(3) == b"BAM"
    except OSError:
        return False


def create_subplots(num_col):
    """
    Prepare the panels for the subplots in case of annotation track.
    """
    if num_col > 1:
        fig = plotly.subplots.make_subplots(
            rows=1,
            cols=num_col,
            column_widths=[0.1, 0.9],
            horizontal_spacing=0.001,
            shared_yaxes=True,
        )
    else:
        fig = plotly.subplots.make_subplots(rows=1, cols=1)
    fig.update_layout(
        {
            "plot_bgcolor": "rgba(0,0,0,0)",
            "paper_bgcolor": "rgba(0,0,0,0)",
        },
        title="methylation frequency",
    )
    return fig


def plot_methylation(subplots, meth_data, num_col):
    samplelist = list(meth_data)
    positionlist = meth_data.index.values.tolist()
    overviewarray = meth_data.to_numpy()

    fig = subplots.add_trace(
        go.Heatmap(z=overviewarray, x=samplelist, y=positionlist), row=1, col=num_col
    )
    fig.update_xaxes(tickangle=45, tickfont=dict(size=10), row=1, col=num_col)
    fig.update_yaxes(tickfont=dict(size=10), row=1, col=num_col)
    return fig


if __name__ == "__main__":
    main()
