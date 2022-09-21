import dash_bio
from scipy.spatial.distance import pdist, squareform
import plotly.figure_factory as ff
import plotly.graph_objects as go
import plotly.graph_objs as go
from plotly.colors import DEFAULT_PLOTLY_COLORS as plcolors
import itertools
from argparse import ArgumentParser
from turtle import title
import pandas as pd
import os
import sys
import numpy as np
import pyranges as pr
from pathlib import Path
import logging
import sys
import subprocess
import gzip
from pathlib import Path
import plotly.express as px
from plotly.subplots import make_subplots


def get_args():
    parser = ArgumentParser(
        description="Create heatmap of methylation frequencies.")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("-f", "--files",
                        nargs="+", help="Haplotype specific methylation frequency files.")
    action.add_argument("-t",
                        "--table", help="Table with methylation frequencies of all samples over region of interest.")
    parser.add_argument("-w", "--window",
                        help="Region to visualise.",
                        required=True if "--example" not in sys.argv else False)  # chr17:44345246-44353106
    parser.add_argument("-n", "--names", nargs="+", default=[])
    parser.add_argument("--gtf", "--gff", help="gtf/gff3 file")
    parser.add_argument(
        "--expand", help="number of base pairs to expand the window with in both directions")
    parser.add_argument("--outtable",
                        help="File to write overview table to.")
    parser.add_argument("--outfig",
                        help="File to write output heatmap to.")
    args = parser.parse_args()
    if len(args.names) == 0:
        for b in args.files:
            c = b.rsplit('.', 2)[0]
            args.names.append(os.path.basename(c))
    if len(args.files) != len(args.names):
        d = len(args.files)
        e = len(args.names)
        sys.exit(
            "Input and names should have same length, length files = %s and length names = %s" % (d, e))
    return args


args = get_args()


class Region(object):
    def __init__(self, region, expand):
        if ':' in region:
            try:
                self.chromosome, interval = region.replace(',', '').split(':')
                try:
                    # see if just integer chromosomes are used
                    self.chromosome = int(self.chromosome)
                except ValueError:
                    pass
                self.begin, self.end = [int(i) for i in interval.split('-')]
                self.begin = self.begin - int(expand)
                if expand:                                  #when expand, make region larger
                    self.end = self.end + int(expand)
                    self.start = self.begin
                self.size = self.end - self.begin
            except ValueError:
                sys.exit("\n\nERROR: Window (-w/--window) inproperly formatted, "
                         "examples of accepted formats are:\n"
                         "'chr5:150200605-150423790' or 'ENST00000647408'\n\n")
            self.size = self.end - self.begin
            if not self.size > 0:
                sys.exit(
                    "\n\nERROR: Window (-w/--window) inproperly formatted, "
                    "begin of the interval has to be smaller than end\n\n"
                )
            self.string = f"{self.chromosome}_{self.begin}_{self.end}"


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

#def main():





#def meth_browser(meth_data, window, expand=False, gft=False, simplify=False):





def get_data(files, names, window):
    """
    Import methylation data from all files in the list methylation_files
    Data can in various formats
    - nanopolish -> calculate_methylation_frequency.py output
    - table in .tsv/.tsv.gz file
    - cram
    data is extracted within the window args.window
    """
    return [read_mods(f, n, window) for f, n in zip(files, names)] #probleem: indien nanopolish_calc_meth_freq input zijn er meerdere files als input, indien overviewtable is er maar 1 file als input


def read_mods(filename, id_name, window, expand):
    """
    converts a file from nanopolish to a pandas dataframe
    input can be from calculate_methylation_frequency
    which will return a dataframe with 'chromosome', 'pos', 'methylated_frequency'
    smoothening the result by a rolling average
    input can also be raw data per read, optionally phased
    which will return a dataframe with 'read', 'chromosome', 'pos', 'log_lik_ratio', 'strand'
    """
    file_type = file_sniffer(filename)
    logging.info(f"File {filename} is of type {file_type}")
    try:
        if file_type.startswith("nanopolish_calc_meth_freq"):
            return parse_nanopolish(filename, id_name, window, expand) 
        elif file_type.startswith("overviewtable")
            return parse_overviewtable()
        # elif file_type in ["cram", "bam"]:
        #     return parse_cram(filename, file_type, name, window) #def parse_cram() moet nog geschreven worden
    except Exception as e:
        logging.error(f"Error processing {filename}.")
        logging.error(e, exc_info=True)
        sys.stderr.write(f"\n\n\nError processing {filename}!\n")
        sys.stderr.write("\n\n\nDetailed error:\n")
        raise


def parse_nanopolish(files, names, window, expand): #input files = calculate_methylation_frequency.py output files
    dfs = []
    region = Region(window, expand)
    for file, name in zip(files, names):
        if window:  #file inlezen met tabix voor gekozen regio
            if Path(file + ".tbi").is_file():
                logging.info(f"Reading {file} using a tabix stream.")
                rg = f"{region.chromosome}:{region.begin}-{region.end}"
                try:
                    tabix_stream = subprocess.Popen(['tabix', file, rg],
                                                    stdout=subprocess.PIPE,
                                                    stderr=subprocess.PIPE)
                except FileNotFoundError as e:
                    logging.error("Error when opening a tabix stream.")
                    logging.error(e, exc_info=True)
                    sys.stderr.write(
                        "\n\nERROR when opening a tabix stream.\n")
                    sys.stderr.write("Is tabix installed and on the PATH?.")
                    raise
                header = gzip.open(file, 'rt').readline().rstrip().split('\t')
                df = pd.read_csv(tabix_stream.stdout,
                                 sep='\t', header=None, names=header)
                logging.info("Read the file in a dataframe.")
        else:
            df = pd.read_csv(file, sep="\t") #telkens hele file inlezen
            logging.info("Read the file in a dataframe.")

        df = pr.PyRanges(df.drop(['group_sequence', 'called_sites_methylated',
                                  'num_motifs_in_group', "called_sites"], axis=1).rename(
            columns={"start": "Start", "chromosome": "Chromosome", "end": "End"}))
        table = df.df
        table["position"] = (table["Start"] + table["End"])/2
        table = table.set_index("position").drop(
            ["Chromosome", "Start", "End"], axis=1)
        return dfs.append(table.rename(columns={'methylated_frequency': name})) #output is an overview table with position as index and for each sample a column with all the methylation frequencies
    
#def parse_overviewtable():
    #input overviewtable inlezen

#def parse_cram():















def args_heatmap(files, window, expand, gtf, simplify):
    dfs = []
    region = Region(window, expand)
    for file, name in zip(files, args.names):
        if window:
            if Path(file + ".tbi").is_file():
                logging.info(f"Reading {file} using a tabix stream.")
                rg = f"{region.chromosome}:{region.begin}-{region.end}"
                try:
                    tabix_stream = subprocess.Popen(['tabix', file, rg],
                                                    stdout=subprocess.PIPE,
                                                    stderr=subprocess.PIPE)
                except FileNotFoundError as e:
                    logging.error("Error when opening a tabix stream.")
                    logging.error(e, exc_info=True)
                    sys.stderr.write(
                        "\n\nERROR when opening a tabix stream.\n")
                    sys.stderr.write("Is tabix installed and on the PATH?.")
                    raise
                header = gzip.open(file, 'rt').readline().rstrip().split('\t')
                df = pd.read_csv(tabix_stream.stdout,
                                 sep='\t', header=None, names=header)
        else:
            df = pd.read_csv(file, sep="\t")

        df = pr.PyRanges(df.drop(['group_sequence', 'called_sites_methylated',
                                  'num_motifs_in_group', "called_sites"], axis=1).rename(
            columns={"start": "Start", "chromosome": "Chromosome", "end": "End"}))
        logging.info("Read the file in a dataframe.")

        table = df.df

        table["MeanPosition"] = (table["Start"] + table["End"])/2

        table = table.set_index("MeanPosition").drop(
            ["Chromosome", "Start", "End"], axis=1)

        dfs.append(table.rename(columns={'methylated_frequency': name}))
    overviewtable = dfs[0].join(dfs[1:], how='outer')
    overviewtable.to_csv(args.outtable, sep="\t", na_rep=np.NaN, header=True)
    samplelist = list(overviewtable)
    positionlist = list(reversed(overviewtable.index.values.tolist()))
    overview = overviewtable.to_numpy()
    overviewarray = np.flip(overview, 0) 
    
    heatmap = make_subplots(rows= 1, cols = 2, column_widths=[0.1, 0.9], horizontal_spacing=0.001, shared_yaxes=True) #indien geen gtf file: optie voor geen subplots toevoegen (moet not "algemener" geschreven worden)
    heatmap.add_trace(go.Heatmap(z=overviewarray, x=samplelist, y=positionlist), row=1, col=2)
    # heatmap.add_trace(go.Heatmap(z=overviewarray, x=samplelist, y=positionlist, text=overviewarray,
    #                   texttemplate="%{text}", textfont={"size": 1}), row=1, col=2) #meth frequencies in de heatmap
    heatmap.update_xaxes(tickangle=45, tickfont=dict(size=10), row=1, col=2) 
    heatmap.update_yaxes(tickfont=dict(size=10),row=1,col=2)

    if gtf:
        annotation_traces, x_pos = gtf_annotation(gtf, region, simplify)
        for annot_trace in annotation_traces:
            heatmap.append_trace(trace=annot_trace, row=1, col=1)
            heatmap.update_xaxes(title_text="", showgrid=False,
                                 showticklabels=False, zeroline=False, row=1, col=1)
            heatmap.update_yaxes(title_text="", showgrid=False,
                                 showticklabels=True, zeroline=False, row=1, col=1)
            heatmap.update_layout({'plot_bgcolor': 'rgba(0,0,0,0)','paper_bgcolor':'rgba(0,0,0,0)',},title="methylation frequency")
    html = heatmap.to_html()  # titel over hele figuur toevoegen

    with open(args.outfig, 'w') as f:
        f.write(html)
    return overviewtable, heatmap

def gtf_annotation(gtf, region, simplify=False):
    result = []
    annotation = parse_annotation(gtf, region, simplify)
    if annotation:
        for x_pos, transcript in enumerate(annotation):
            line = make_per_gene_annot_line_trace(transcript, region, x_pos)
            exons = [make_per_exon_arrow_trace(transcript, begin, end, x_pos)
                     for begin, end in transcript.exon_tuples
                     if region.begin < begin and region.end > end]
            result.extend([line, *exons])
        return result, x_pos
    else:
        return result, 0


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
        df["begin"].between(region.begin, region.end) | df["end"].between(
            region.begin, region.end),
        feature,
    ].unique()


def assign_colors_to_genes(transcripts):
    genes = set([t.gene for t in transcripts])
    colordict = {g: c for g, c in zip(genes, plcolors * 100)}
    for t in transcripts:
        t.color = colordict[t.gene]


def parse_annotation(gff, region, simplify=False):
    type = annot_file_type(gff)
    logging.info(f"Parsing {type} file...")
    """
    Parse the gff and select the relevant region as determined by the window
    return as Transcript objects
    """
    if Path(gff + ".tbi").is_file():
        logging.info(f"Reading {gff} using a tabix stream.")
        rg = f"{region.chromosome}:{region.begin}-{region.end}"
        try:
            tabix_stream = subprocess.Popen(['tabix', gff, rg],
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE)
        except FileNotFoundError as e:
            logging.error("Error when opening a tabix stream.")
            logging.error(e, exc_info=True)
            sys.stderr.write(
                "\n\nERROR when opening a tabix stream.\n")
            sys.stderr.write("Is tabix installed and on the PATH?.")
            raise
    annotationfile = pd.read_csv(tabix_stream.stdout, sep='\t', header=None, names=[
                                 "chromosome", "source", "feature_type", "begin", "end", "score", "strand", "phase", "attributes"])
    #else: #######

    annotationfile.to_csv(
        "/home/ecoopman/outputresults/annotationtabletest.tsv", sep="\t", na_rep=np.NaN, header=True)

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
            transcript.append('')

    annotationfile["gene"] = gene
    annotationfile["transcript"] = transcript

    annotationfile.drop(["source", "feature_type", "score",
                        "phase", "attributes"], axis=1, inplace=True)

    if simplify:
        annotationfile.drop_duplicates(
            subset=["chromosome", "begin", "end", "gene"], inplace=True)
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
    return go.Scatter(y=[max(transcript.begin, region.begin),
                         min(transcript.end, region.end)],
                      x=[x_pos, x_pos],
                      mode='lines',
                      line=dict(width=2, color=transcript.color),
                      name=transcript.transcript,
                      text=transcript.gene,
                      hoverinfo='text',
                      showlegend=False)


def make_per_exon_arrow_trace(transcript, begin, end, x_pos):
    """Generate a line+marker trace for the exon
    The shape is an arrow, as defined by the strand in transcript.marker
    """
    return go.Scatter(y=[begin, end],
                      x=[x_pos, x_pos],
                      mode='lines+markers',
                      line=dict(width=8, color=transcript.color),
                      name=transcript.transcript,
                      text=transcript.gene,
                      hoverinfo='text',
                      showlegend=False,
                      marker=dict(symbol=transcript.marker,
                                  size=8))


def file_sniffer(filename):
    """
    Takes in a filename and tries to guess the input file type
    """
    if not Path(filename).is_file():
        sys.exit(
            f"\n\nERROR: File {filename} does not exist, please check the path!\n")
    if is_bam_file(filename): #input BAM
        return "bam"
    if is_cram_file(filename): #input CRAM
        return "cram"
    if is_gz_file(filename): #input: calculate_methylation_frequency.py output (.tsv or .tsv.gz) OR own overview table (.tsv or .tsv.gz)
        import gzip

        header = gzip.open(filename, "rt").readline()
    else:
        header = open(filename, "r").readline()

    if "num_motifs_in_group" in header:
        return "nanopolish_calc_meth_freq" #calculate_methylation_frequency.py output as input
    if "position" in header:
        return "overviewtable" #own overviewtable as input: needs first column header to be "position"
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

overviewtable, heatmap = (
    args.files, args.names, args.window, args.expand, args.gtf, simplify=False)
