import itertools
import logging
import shlex
import sys
import pandas as pd
from pathlib import Path
from methylmap.region import Region
import subprocess
from plotly.colors import DEFAULT_PLOTLY_COLORS as plcolors
import plotly.graph_objects as go


class Transcript(object):
    def __init__(self, transcript, gene, exon_tuples, strand):
        self.transcript = transcript
        self.gene = gene
        self.exon_tuples = list(exon_tuples)
        self.strand = strand
        self.marker = (
            "triangle-down" if self.strand == "+" else "triangle-up"
        )  # reversed
        self.begin = min(list(itertools.chain.from_iterable(self.exon_tuples)))
        self.end = max(list(itertools.chain.from_iterable(self.exon_tuples)))
        self.color = ""


def annot_file_type(annot_file):
    """
    Figure out type of annotation file.
    """
    if annot_file.endswith(
        (".gff.gz",".gff2.gz", ".gff3.gz")
    ):
        return "gff"
    elif annot_file.endswith((".gff",".gff2",".gff3")):
        logging.error("Annotation file not bgzipped.")
        sys.exit(
            "ERROR: annotation file not bgzipped.\n"
        )
    else:
        logging.error("Unrecognized extension of the annotation file. Supported are gff.gz, gff2.gz, and gff3.gz")
        sys.exit(
            "ERROR: unrecognized extension of the annotation file.\n"
            "Supported are gff.gz, gff2.gz, and gff3.gz"
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


def parse_annotation(gff, window):
    """
    Parse the gff and select the relevant window as determined by the window input by using tabix
    """
    type = annot_file_type(gff)
    logging.info(f"Parsing {type} file...")
    if not Path(gff + ".tbi").is_file():
        try:
            logging.info(
                "Make .tbi file from annotion file for fast selection window of interest."
            )
            tabix_gff = subprocess.Popen(
                shlex.split(f"tabix -p gff {gff}"), stderr=subprocess.PIPE
            )
        except FileNotFoundError as e:
            logging.error("Error when making a .tbi file.")
            logging.error(e, exc_info=True)
            sys.stderr.write("\n\nERROR when making a .tbi file.\n")
            sys.stderr.write("Is tabix installed and on the PATH?.")
            sys.stderr.write(f"\n\n\nDetailed error: {tabix_gff.stderr.read()}\n")
            raise
        if tabix_gff.returncode:
            sys.exit(f"\n\n\nReceived tabix error:\n{tabix_gff.stderr.read()}\n")
    try:
        logging.info(f"Reading {gff} using a tabix stream.")
        tabix_stream = subprocess.Popen(
            ["tabix", gff, window.fmt], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    except FileNotFoundError as e:
        logging.error("Error when opening a tabix stream.")
        logging.error(e, exc_info=True)
        sys.stderr.write("\n\nERROR when opening a tabix stream.\n")
        sys.stderr.write(f"\n\n\nDetailed error: {tabix_stream.stderr.read()}\n")
        raise
    if tabix_stream.returncode:
        sys.exit(f"\n\n\nReceived tabix error:\n{tabix_stream.stderr.read()}\n")
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
        sys.stderr.write(f"\nFound {len(res)} gene(s) in the window.\n")
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


def make_per_gene_annot_line_trace(transcript, window, x_pos):
    """Generate a line trace for the gene.
    Trace can get limited by the window sizes.
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
