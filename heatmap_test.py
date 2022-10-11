import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots
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
        self.start = (
            self.begin
        )  # start is now just an alias for begin because I tend to forget
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
        cases=args.cases,
        controls=args.controls,
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
    parser.add_argument("-n", "--names", nargs="+", default=[])
    parser.add_argument("--gff", "--gtf", help="gtf/gff3 file")
    parser.add_argument(
        "--expand",
        help="number of base pairs to expand the window with in both directions",
        type=int,
        default=0,
    )
    parser.add_argument("--outtable", help="File to write the frequencies table to.")
    parser.add_argument("--outfig", help="File to write output heatmap to.")
    parser.add_argument("--groups",help="List of category per sample")
    parser.add_argument("-s", "--simplify", action="store_true")  # default: False
    args = parser.parse_args()
    # werkt niet wanneer frequencies table als input: args.names input kan niet van basename files gehaald worden; ofwel list met names als input, ofwel optie om de header als names te gebruiken, moet nog geschreven worden
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
    if window:  # and expand?
        window = Region(window, expand)
    heatmapfig = 1
    if gff:
        annotationfig = 1
    else:
        annotationfig = 0
    num_col = heatmapfig + annotationfig  # number of subplots (columns) needed
    subplots = create_subplots(num_col)

    # frequencies table with all meth frequencies of all samples
    meth_data, window = read_mods(files, table, names, window, groups, gff, outtable)
    fig = plot_methylation(subplots, meth_data, num_col)

    if gff:
        annotation_traces = gff_annotation(gff, window, simplify)
        for annot_trace in annotation_traces:
            fig.append_trace(trace=annot_trace, row=1, col=1)
        fig.update_xaxes(
            title_text="", showticklabels=False, zeroline=False, row=1, col=1
        )
        fig.update_yaxes(
            title_text="", showticklabels=True, zeroline=False, row=1, col=1
        )
    with open(outfig, "w") as f:
        f.write(fig.to_html())
    return meth_data, fig


def read_mods(files, table, names, window, groups, gff, outtable):
    """
    converts a file from nanopolish to a pandas dataframe
    input can be from calculate_methylation_frequency
    which will return a dataframe with 'chromosome', 'pos', 'methylated_frequency'
    """
    # checkt alleen eerste file in de list met alle files, niet ideaal, hoogstwss ook probleem wanneer enkel frequencies table als input o
    # checkt alleen eerste file in de list met alle files
    if files:
        file_type = file_sniffer(files[0])
    elif table:
        file_type = file_sniffer(table)
    # logging.info(f"File {files[0]} is of type {file_type}")
    try:
        if file_type == "nanopolish_calc_meth_freq":
            return parse_nanopolish(files, table, names, window, groups, outtable)
        elif file_type == "methfrequencytable": 
            return parse_methfrequencytable(table, names, window, groups, gff, outtable)
        # elif file_type == "overviewtable":
        #     return parse_overviewtable(table, outtable)
        # elif file_type in ["cram", "bam"]: ####def parse_cram() moet nog geschreven worden
        #     return parse_cram(filename, file_type, name, window)
    except Exception as e:
        # WDC: I don't think it is necessarily correct that files[0] is the cause
        # logging.error(f"Error processing {files[0]}.")
        # logging.error(e, exc_info=True)
        # sys.stderr.write(f"\n\n\nError processing {files[0]}!\n")
        sys.stderr.write("\n\n\nDetailed error:\n")
        raise


def parse_nanopolish(files, table, names, window, groups, outtable):
    """
    input files = calculate_methylation_frequency.py output files
    """
    if table:
        overviewtable = pd.read_table(table).sort_values(["group","name"])
        files = overviewtable["path"].tolist()
        names = overviewtable["name"].tolist()

    dfs = []
    for file, name in zip(files, names):
        if window:  # file inlezen met tabix voor gekozen regio
            if Path(
                file + ".tbi"
            ).is_file():  # nog optie schrijven om tbi file te maken
                logging.info(f"Reading {file} using a tabix stream.")
                try:
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
                table = pd.read_csv(
                    tabix_stream.stdout, sep="\t", header=None, names=header
                )
                logging.info("Read the file in a dataframe.")
            # else: # zelf .tbi file maken #####Dit stuk code werkt nog niet: hoe teruggaan naar 'if' hierboven wanneer deze else is uitgevoerd?
            #     try:
            #         with open(file) as f:
            #             subprocess.Popen(["tabix","-S1", "-s1", "-b2", "-e3", file], stdout=f)
            #     except FileNotFoundError as e:
            #         logging.error("Error when making a .tbi file.")
            #         logging.error(e, exc_info=True)
            #         sys.stderr.write(
            #             "\n\nERROR when making a .tbi file.\n")
            #         sys.stderr.write("Is tabix installed and on the PATH?.")
            #         raise
            #     #return parse_nanopolish(files, names, window) #geen goed idee, telkens volledig opnieuw beginnen, beter naar 'if' teruggaan dan helemaal naar begin
            # WDC change the logic and check if the tbi does not exist, then make it, and continue with the code for reading from tabix_stream

        # probleem wanneer én geen .tbi file beschikbaar én van file kan geen .tbi file gemaakt worden
        else:
            table = pd.read_csv(file, sep="\t")  # telkens hele file inlezen
            logging.info("Read the file in a dataframe.")
            #if gff: also extract window out of it? 

        df = pr.PyRanges(
            table.drop(
                [
                    "group_sequence",
                    "called_sites_methylated",
                    "num_motifs_in_group",
                    "called_sites",
                ],
                axis=1,
            ).rename(
                columns={"start": "Start", "chromosome": "Chromosome", "end": "End"}
            )
        )
        table = df.df
        table["position"] = (table["Start"] + table["End"]) / 2
        table = table.set_index("position").drop(["Chromosome", "Start", "End"], axis=1)
        dfs.append(table.rename(columns={"methylated_frequency": name}))
    methfrequencytable = dfs[0].join(dfs[1:], how="outer")
    reversedtable = methfrequencytable.reindex(
        index=methfrequencytable.index[::-1]
    ) 
    # output is an meth frequency table with position as index and for each sample a column with all the methylation frequencies
    
    if groups and not table: #kan dit zo?
        headerlist = list(reversedtable.columns.values)
        res = []
        for x in zip(headerlist, groups):
            res.append(x)
            output = sorted(res, key=lambda x: x[-1])
            orderedlist = []
            for i in output:
                orderedlist.append(i[0])
            methfreqtable = reversedtable.reindex(columns=orderedlist)

    methfreqtable.to_csv(outtable, sep="\t", na_rep=np.NaN, header=True)
    return methfreqtable


def parse_methfrequencytable(
    table, names, window, groups, gff, outtable
):
    # if expand: #kan niet in args want --table kan ook overviewtable als input krijgen, waar --extend wel mogelijk is
    #         logging.warning("Error when reading in methylation frequency table: no --extend option possible.")
    #         sys.stderr.write("\n\nERROR when reading in methylation frequency tabloe:\n")
    #         sys.stderr.write("--extend option not possible.")
    table = pd.read_table(table)
    for x in table["position"]:
        if ":" in x:
            table[['chrom', 'position']] = table['position'].str.split(':', 1, expand=True)
    
    numberofpositions = len(table) - 1
    value1 = table.iloc[0, 0]
    value2 = table.iloc[numberofpositions, 0]
    if (
        value2 > value1
    ):  # reverse table when last position is further than first position #####wat met regio's die over verschillende chromosomen gaan???
        methfreqtable = table.reindex(index=table.index[::-1])
    else:
        methfreqtable = table
    
    if window:  # select window out of table
        methfreqtable = methfreqtable[
            (methfreqtable["position"].astype(float) >= window.begin)
            & (methfreqtable["position"].astype(float) <= window.end)
        ]
    else:
        if gff: #if annotation and no --window option present, extract window out of methfreqtable
            if "chrom" in methfreqtable.columns:
                begin = float(methfreqtable.iloc[numberofpositions,0])
                end = float(methfreqtable.iloc[0,0])
                if (methfreqtable['chrom'] == methfreqtable['chrom'][0]).all():
                    chrom = methfreqtable.iloc[0, methfreqtable.columns.get_loc('chrom')]
                    if not chrom.startswith("chr"):
                        logging.error("Error when extracting window out of methfreqtable positions.")
                        sys.stderr.write("\n\nError when extracting window out of methfreqtable positions.\n")
                        sys.stderr.write("Is position column of methfreqtable in format chr1:123456?") 
                else:
                    logging.error("Error when extracting window out of methfreqtable positions.")
                    sys.stderr.write("\n\nError when extracting window out of methfreqtable positions.\n")
                    sys.stderr.write("Window over different chromosomes is not possible.")
            window = Region(f"{chrom}:{round(begin)}-{round(end)}")

    methfreqtable.drop(["chrom"], axis=1, inplace=True)
    methfreqtable.set_index("position", inplace=True)
    
    if names:
        methfreqtable.columns = names
        sys.stderr.write("Error when matching the input names to the input table")
        sys.stderr.write(
            "Is the number of the input names same as the number of samples in the table?."
        )
        logging.info("Error when matching the input names to the input table")

    if groups:
        try: 
            headerlist = list(methfreqtable.columns.values)
            res = []
            for x in zip(headerlist, groups):
                res.append(x)
                output = sorted(res, key=lambda x: x[-1])
                orderedlist = []
                for i in output:
                    orderedlist.append(i[0])
                methfreqtable = methfreqtable.reindex(columns=orderedlist)
        except FileNotFoundError as e:
            logging.error("Error when matching --groups with samples.")
            logging.error(e, exc_info=True)
            sys.stderr.write("\n\nError when matching --groups with samples.\n")
            sys.stderr.write(f"Is length of the --groups list ({len(groups)}) matching with number of sample files?")
            raise
    
    methfreqtable.to_csv(outtable, sep="\t", na_rep=np.NaN, header=True)
    return methfreqtable, window


# def parse_cram():
# input cram files inlezen


def gff_annotation(gff, window, simplify):
    result = []
    annotation = parse_annotation(gff, window, simplify)
    # Isn't this always True? The gff_annotation function is only called if args.gff is not None
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
    elif annot_file.endswith(
        (".gff", ".gff.gz", ".gff2", ".gff2.gz", ".gff3", ".gff3.gz")
    ):
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
        df["begin"].between(window.begin, window.end)
        | df["end"].between(window.begin, window.end),
        feature,
    ].unique()


def assign_colors_to_genes(transcripts):
    genes = set([t.gene for t in transcripts])
    colordict = {g: c for g, c in zip(genes, plcolors * 100)}
    for t in transcripts:
        t.color = colordict[t.gene]


# WDC This function is too long
def parse_annotation(gff, window, simplify):
    """
    Parse the gff and select the relevant window as determined by the window
    return as Transcript objects
    """
    type = annot_file_type(gff)
    logging.info(f"Parsing {type} file...")
    if Path(gff + ".tbi").is_file():  # nog optie schrijven om tbi file te maken
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
    # zelf .tbi file maken en terug naar begin (parse_annotation()) gaan, dit stuk code werkt!
    else:
        try:
            with open(gff) as f:
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
    # probleem wanneer én geen .tbi file beschikbaar én van file kan geen .tbi file gemaakt worden; annotationfile zonder tabix inlezen en dan window selecteren?

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
        annotationfile.drop_duplicates(
            subset=["chromosome", "begin", "end", "gene"], inplace=True
        )
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
    if "path" in header:
        df = pd.read(filename)
        if df['path'].iloc[0].endswith(".tsv"): #only checks first file in overviewtable
            return "nanopolish_calc_meth_freq"
        elif df['path'].iloc[0].endswith(".tsv.gz"): #only checks first file in overviewtable
            return "nanoplish_calc_meth_freq"
        elif df['path'].iloc[0].endswith(".bam"): #only checks first file in overviewtable
            return "bam"
        elif df['path'].iloc[0].endswith(".cram"): #only checks first file in overviewtable
            return "cram"
    if header.startswith("position"):
        # own methfrequencytable as input: needs first column header to be "position" (nog niet de ideale oplossing)
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
