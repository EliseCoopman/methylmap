import logging
import pandas as pd
from pathlib import Path
import subprocess
from plotly.colors import DEFAULT_PLOTLY_COLORS as plcolors



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
    if not Path(gff + ".tbi").is_file():
        try:
            with open(gff):
                ### WDC Does this attempt to overwrite the gff?
                ### WDC is there even a stdout with this command?
                subprocess.Popen(
                    ["tabix", "-s1", "-b4", "-e5", gff],stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                ) 
                # possible problem: -S, -s, -b, -e zijn anders in andere gff files
                # WDC "tabix ,"-p", gff" -> dit werkt niet
        except FileNotFoundError as e:
            logging.error("Error when making a .tbi file.")
            logging.error(e, exc_info=True)
            sys.stderr.write("\n\nERROR when making a .tbi file.\n")
            sys.stderr.write("Is tabix installed and on the PATH?.")
            raise     
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
        raise
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

class Transcript(object):
    def __init__(self, transcript, gene, exon_tuples, strand):
        self.transcript = transcript
        self.gene = gene
        self.exon_tuples = list(exon_tuples)
        self.strand = strand
        self.marker = "triangle-down" if self.strand == "+" else "triangle-up" #reversed
        self.begin = min(list(itertools.chain.from_iterable(self.exon_tuples)))
        self.end = max(list(itertools.chain.from_iterable(self.exon_tuples)))
        self.color = ""