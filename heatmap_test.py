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
    def __init__(self, region, expand):
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
        self.fmt = f"{region.chromosome}:{region.begin}-{region.end}"


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
    # WDC returned meth_data and fig are not used, probably not necessary to return those?
    meth_data, fig = meth_browser(
        files=args.files,
        outtable=args.outtable,
        outfig=args.outfig,
        names=args.names,
        window=args.window,
        expand=args.expand,
        gff=args.gff,
        simplify=args.simplify,
    )


def get_args():
    parser = ArgumentParser(description="Create heatmap of methylation frequencies.")
    parser.add_argument(
        "-f", "--files", nargs="+", help="Haplotype specific methylation frequency files."
    )
    # WDC you don't have a --example (yet)
    # chr17:44345246-44353106
    parser.add_argument(
        "-w",
        "--window",
        help="Region to visualise.",
        required=True if "--example" not in sys.argv else False,
    )
    parser.add_argument("-n", "--names", nargs="+", default=[])
    parser.add_argument("--gff", "--gtf", help="gtf/gff3 file")
    parser.add_argument(
        "--expand",
        help="number of base pairs to expand the window with in both directions",
        type=int,
        default=0,
    )
    parser.add_argument("--outtable", help="File to write overview table to.")
    parser.add_argument("--outfig", help="File to write output heatmap to.")
    args = parser.parse_args()
    # werkt niet wanneer overview table als input: args.names input kan niet van basename files gehaald worden; ofwel list met names als input, ofwel optie om de header als names te gebruiken, moet nog geschreven worden
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
    outfig="browser.html",
    outtable=False,
    names=False,
    window=False,
    expand=False,
    gff=False,
    simplify=False,
):
    if window:
        region = Region(window, expand)
    heatmapfig = 1
    if gff:
        annotationfig = 1
    else:
        annotationfig = 0
    num_col = heatmapfig + annotationfig  # number of subplots (columns) needed
    subplots = create_subplots(num_col)

    # overviewtable with all meth frequencies of all samples
    meth_data = read_mods(files, names, region, outtable)
    fig = plot_methylation(subplots, meth_data, num_col)

    if gff:
        annotation_traces = gff_annotation(gff, region, simplify)
        for annot_trace in annotation_traces:
            fig.append_trace(trace=annot_trace, row=1, col=1)
        fig.update_xaxes(title_text="", showticklabels=False, zeroline=False, row=1, col=1)
        fig.update_yaxes(title_text="", showticklabels=True, zeroline=False, row=1, col=1)
    with open(outfig, "w") as f:
        f.write(fig.to_html())
    return meth_data, fig


def read_mods(files, names, region, outtable):
    """
    converts a file from nanopolish to a pandas dataframe
    input can be from calculate_methylation_frequency
    which will return a dataframe with 'chromosome', 'pos', 'methylated_frequency'
    """
    # checkt alleen eerste file in de list met alle files, niet ideaal, hoogstwss ook probleem wanneer enkel overviewtable als input o
    file_type = file_sniffer(files[0])  # checkt alleen eerste file in de list met alle files
    logging.info(f"File {files[0]} is of type {file_type}")
    try:
        if file_type == "nanopolish_calc_meth_freq":
            return parse_nanopolish(files, names, region, outtable)
        # elif file_type == "overviewtable":  #nog niet volledig werkende
        #     return parse_overviewtable(files)
        # elif file_type in ["cram", "bam"]: ####def parse_cram() moet nog geschreven worden
        #     return parse_cram(filename, file_type, name, window)
    except Exception as e:
        # WDC: I don't think it is necessarily correct that files[0] is the cause
        logging.error(f"Error processing {files[0]}.")
        logging.error(e, exc_info=True)
        sys.stderr.write(f"\n\n\nError processing {files[0]}!\n")
        sys.stderr.write("\n\n\nDetailed error:\n")
        raise


def parse_nanopolish(files, names, region, outtable):
    """
    input files = calculate_methylation_frequency.py output files
    """
    dfs = []
    for file, name in zip(files, names):
        if region:  # file inlezen met tabix voor gekozen regio
            if Path(file + ".tbi").is_file():  # nog optie schrijven om tbi file te maken
                logging.info(f"Reading {file} using a tabix stream.")
                try:
                    tabix_stream = subprocess.Popen(
                        ["tabix", file, region.fmt], stdout=subprocess.PIPE, stderr=subprocess.PIPE
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
            #     #return parse_nanopolish(files, names, region) #geen goed idee, telkens volledig opnieuw beginnen, beter naar 'if' teruggaan dan helemaal naar begin
            # WDC change the logic and check if the tbi does not exist, then make it, and continue with the code for reading from tabix_stream

        # probleem wanneer én geen .tbi file beschikbaar én van file kan geen .tbi file gemaakt worden
        else:
            table = pd.read_csv(file, sep="\t")  # telkens hele file inlezen
            logging.info("Read the file in a dataframe.")

        df = pr.PyRanges(
            table.drop(
                [
                    "group_sequence",
                    "called_sites_methylated",
                    "num_motifs_in_group",
                    "called_sites",
                ],
                axis=1,
            ).rename(columns={"start": "Start", "chromosome": "Chromosome", "end": "End"})
        )
        table = df.df
        table["position"] = (table["Start"] + table["End"]) / 2
        table = table.set_index("position").drop(["Chromosome", "Start", "End"], axis=1)
        dfs.append(table.rename(columns={"methylated_frequency": name}))
    overviewtable = dfs[0].join(dfs[1:], how="outer")
    # output is an overview table with position as index and for each sample a column with all the methylation frequencies
    overviewtable.to_csv(outtable, sep="\t", na_rep=np.NaN, header=True)
    return overviewtable


# def parse_overviewtable(files): #overviewtable input option not yet working
#     return pd.read_table(files)

# def parse_cram():
# input cram files inlezen


def gff_annotation(gff, region, simplify):
    result = []
    annotation = parse_annotation(gff, region, simplify)
    # Isn't this always True? The gff_annotation function is only called if args.gff is not None
    if annotation:
        for x_pos, transcript in enumerate(annotation):
            line = make_per_gene_annot_line_trace(transcript, region, x_pos)
            exons = [
                make_per_exon_arrow_trace(transcript, begin, end, x_pos)
                for begin, end in transcript.exon_tuples
                if region.begin < begin and region.end > end
            ]
            result.extend([line, *exons])
        return result
    else:
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


def transcripts_in_window(df, region, feature="transcript"):
    """
    Return the transcript names for which
    either the end or the begin of an exon is within the window
    """
    return df.loc[
        df["begin"].between(region.begin, region.end) | df["end"].between(region.begin, region.end),
        feature,
    ].unique()


def assign_colors_to_genes(transcripts):
    genes = set([t.gene for t in transcripts])
    colordict = {g: c for g, c in zip(genes, plcolors * 100)}
    for t in transcripts:
        t.color = colordict[t.gene]


# WDC This function is too long
def parse_annotation(gff, region, simplify):
    """
    Parse the gff and select the relevant region as determined by the window
    return as Transcript objects
    """
    type = annot_file_type(gff)
    logging.info(f"Parsing {type} file...")
    if Path(gff + ".tbi").is_file():  # nog optie schrijven om tbi file te maken
        logging.info(f"Reading {gff} using a tabix stream.")
        try:
            tabix_stream = subprocess.Popen(
                ["tabix", gff, region.fmt], stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
        except FileNotFoundError as e:
            logging.error("Error when opening a tabix stream.")
            logging.error(e, exc_info=True)
            sys.stderr.write("\n\nERROR when opening a tabix stream.\n")
            sys.stderr.write("Is tabix installed and on the PATH?.")
            raise
    else:  # zelf .tbi file maken en terug naar begin (parse_annotation()) gaan, dit stuk code werkt!
        try:
            with open(gff) as f:
                subprocess.Popen(
                    ["tabix", "-s1", "-b4", "-e5", gff], stdout=f
                )  # possible problem: -S, -s, -b, -e zijn anders in andere gff files
                # WDC "tabix -p gff"
        except FileNotFoundError as e:
            logging.error("Error when making a .tbi file.")
            logging.error(e, exc_info=True)
            sys.stderr.write("\n\nERROR when making a .tbi file.\n")
            sys.stderr.write("Is tabix installed and on the PATH?.")
            raise
        return parse_annotation(gff, region, simplify=False)
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
        annotationfile.drop_duplicates(subset=["chromosome", "begin", "end", "gene"], inplace=True)
        res = []
        for g in transcripts_in_window(annotationfile, region, feature="gene"):
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
        sys.stderr.write(f"Found {len(res)} gene(s) in the region.\n")
        logging.info(f"Found {len(res)} gene(s) in the region.\n")
    else:
        res = []
        for t in transcripts_in_window(annotationfile, region, feature="transcript"):
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
        sys.stderr.write(f"Found {len(res)} transcript(s) in the region.\n")
        logging.info(f"Found {len(res)} transcript(s) in the region.\n")
    assign_colors_to_genes(res)
    return res


def make_per_gene_annot_line_trace(transcript, region, x_pos):
    """Generate a line trace for the gene
    Trace can get limited by the window sizes
    """
    return go.Scatter(
        y=[max(transcript.begin, region.begin), min(transcript.end, region.end)],
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
    # input: calculate_methylation_frequency.py output (.tsv or .tsv.gz) OR own overview table (.tsv or .tsv.gz)
    if is_gz_file(filename):
        import gzip

        header = gzip.open(filename, "rt").readline()
    else:
        header = open(filename, "r").readline()

    if "num_motifs_in_group" in header:
        # calculate_methylation_frequency.py output as input
        return "nanopolish_calc_meth_freq"
    if header.startswith("position"):
        # own overviewtable as input: needs first column header to be "position" (nog niet de ideale oplossing)
        return "overviewtable"
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
    positionlist = list(reversed(meth_data.index.values.tolist()))  # WDC reversed?
    overview = meth_data.to_numpy()
    overviewarray = np.flip(overview, 0)

    fig = subplots.add_trace(
        go.Heatmap(z=overviewarray, x=samplelist, y=positionlist), row=1, col=num_col
    )
    fig.update_xaxes(tickangle=45, tickfont=dict(size=10), row=1, col=num_col)
    fig.update_yaxes(tickfont=dict(size=10), row=1, col=num_col)
    return fig


if __name__ == "__main__":
    main()
