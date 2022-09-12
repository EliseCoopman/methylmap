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
    def __init__(self, region, fasta=None):
        if ':' in region:
            try:
                self.chromosome, interval = region.replace(',', '').split(':')
                try:
                    # see if just integer chromosomes are used
                    self.chromosome = int(self.chromosome)
                except ValueError:
                    pass
                self.begin, self.end = [int(i) for i in interval.split('-')]
                self.start = self.begin
                self.size = self.end - self.begin
                # if self.size < 100000:
                #     self.begin = self.begin - self.size
                #     self.end = self.end + self.size
                # else:
                #     self.begin = self.begin - 100000
                #     self.end = self.end + 100000
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


def args_heatmap(files, window):
    dfs = []
    region = Region(window)
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
    joined = dfs[0].join(dfs[1:], how='outer')
    joined.to_csv(args.outtable, sep="\t", na_rep=np.NaN, header=True)
    heatmap = px.imshow(joined, text_auto=True)
    heatmap.update_xaxes(tickangle=45, tickfont=dict(size=8))
    heatmap.update_yaxes(tickfont=dict(size=10))
    heatmap.update_layout(
        title="#######")
    html = heatmap.to_html()

    with open(args.outfig, 'w') as f:
        f.write(html)
    return joined, heatmap


joined, heatmap = args_heatmap(args.files, args.window)

###### heatmap-ANNOTATION

# def get_args():
#     parser = ArgumentParser(
#         description="Create heatmap of methylation frequency")
#     action = parser.add_mutually_exclusive_group(required=True)
#     action.add_argument("--files",
#                         nargs="+", help="Haplotype specific frequency files.")
#     action.add_argument(
#         "--joinedtable", help="Table with combined information")
#     parser.add_argument("-w", "--window",
#                         help="window (region) to which the visualisation has to be restricted",
#                         required=True if "--example" not in sys.argv else False)  # chr17:44345246-44353106
#     parser.add_argument("--gff", help="gff3 file")
#     parser.add_argument("--names", nargs="+", default=[])
#     parser.add_argument("--outtable",
#                         help="File to write results to.")
#     parser.add_argument("--outfig",
#                         help="File to write results to.")
#     args = parser.parse_args()
#     if len(args.names) == 0:
#         for b in args.files:
#             c = b.rsplit('.', 2)[0]
#             args.names.append(os.path.basename(c))
#     if len(args.files) != len(args.names):
#         d = len(args.files)
#         e = len(args.names)
#         sys.exit(
#             "input and names should have same length, length files = %s and length names = %s" % (d, e))
#     return args

# args = get_args()

# class Region(object):
#     def __init__(self, region, fasta=None):
#         if ':' in region:
#             try:
#                 self.chromosome, interval = region.replace(',', '').split(':')
#                 try:
#                     # see if just integer chromosomes are used
#                     self.chromosome = int(self.chromosome)
#                 except ValueError:
#                     pass
#                 self.begin, self.end = [int(i) for i in interval.split('-')]
#                 self.start = self.begin
#             except ValueError:
#                 sys.exit("\n\nERROR: Window (-w/--window) inproperly formatted, "
#                          "examples of accepted formats are:\n"
#                          "'chr5:150200605-150423790' or 'ENST00000647408'\n\n")
#             self.size = self.end - self.begin
#             if not self.size > 0:
#                 sys.exit(
#                     "\n\nERROR: Window (-w/--window) inproperly formatted, "
#                     "begin of the interval has to be smaller than end\n\n"
#                 )
#             self.string = f"{self.chromosome}_{self.begin}_{self.end}"


# class Transcript(object):
#     def __init__(self, transcript, gene, exon_tuples, strand):
#         self.transcript = transcript
#         self.gene = gene
#         self.exon_tuples = list(exon_tuples)
#         self.strand = strand
#         self.marker = "triangle-right" if self.strand == "+" else "triangle-left"
#         self.begin = min(list(itertools.chain.from_iterable(self.exon_tuples)))
#         self.end = max(list(itertools.chain.from_iterable(self.exon_tuples)))
#         self.color = ""


# def args_heatmap(files, window, gtf, simplify):
#     dfs = []
#     region = Region(window)
#     for file, name in zip(files, args.names):
#         if window:
#             if Path(file + ".tbi").is_file():
#                 logging.info(f"Reading {file} using a tabix stream.")
#                 rg = f"{region.chromosome}:{region.begin}-{region.end}"
#                 try:
#                     tabix_stream = subprocess.Popen(['tabix', file, rg],
#                                                     stdout=subprocess.PIPE,
#                                                     stderr=subprocess.PIPE)
#                 except FileNotFoundError as e:
#                     logging.error("Error when opening a tabix stream.")
#                     logging.error(e, exc_info=True)
#                     sys.stderr.write(
#                         "\n\nERROR when opening a tabix stream.\n")
#                     sys.stderr.write("Is tabix installed and on the PATH?.")
#                     raise
#                 header = gzip.open(file, 'rt').readline().rstrip().split('\t')
#                 df = pd.read_csv(tabix_stream.stdout,
#                                  sep='\t', header=None, names=header)
#         else:
#             df = pd.read_csv(file, sep="\t")

#         df = pr.PyRanges(df.drop(['group_sequence', 'called_sites_methylated',
#                                   'num_motifs_in_group', "called_sites"], axis=1).rename(
#             columns={"start": "Start", "chromosome": "Chromosome", "end": "End"}))
#         logging.info("Read the file in a dataframe.")

#         table = df.df

#         table["MeanPosition"] = (table["Start"] + table["End"])/2

#         table = table.set_index("MeanPosition").drop(
#             ["Chromosome", "Start", "End"], axis=1)

#         dfs.append(table.rename(columns={'methylated_frequency': name}))
#     joined = dfs[0].join(dfs[1:], how='outer')

#     joined.to_csv(args.outtable, sep="\t", na_rep=np.NaN, header=True)

#     heatmap = px.imshow(joined, text_auto=True)
#     heatmap.update_xaxes(tickangle=45, tickfont=dict(size=10))
#     heatmap.update_yaxes(tickfont=dict(size=10))

#     if gtf:
#         sys.stderr.write("0")
#         annotation_traces, y_max = gtf_annotation(gtf, region, simplify)
#         for annot_trace in annotation_traces:
#             heatmap.append_trace(trace=annot_trace, row=1, col=1)

#     html = heatmap.to_html()

#     with open(args.outfig, 'w') as f:
#         f.write(html)
#     return joined, heatmap


# def gtf_annotation(gtf, region, simplify=False):
#     result = []
#     annotation = parse_annotation(gtf, region, simplify)
#     if annotation:
#         for y_pos, transcript in enumerate(annotation):
#             line = make_per_gene_annot_line_trace(transcript, region, y_pos)
#             exons = [make_per_exon_arrow_trace(transcript, begin, end, y_pos)
#                      for begin, end in transcript.exon_tuples
#                      if region.begin < begin and region.end > end]
#             result.extend([line, *exons])
#         return result, y_pos
#     else:
#         return result, 0


# def annot_file_sniffer(annot_file):
#     """
#     Figure out type of annotation file
#     Right not just lazily focus on the extension
#     """
#     if annot_file.endswith((".gtf", ".gtf.gz")):
#         return "gtf"
#     elif annot_file.endswith((".gff", ".gff.gz", ".gff2", ".gff2.gz", ".gff3", ".gff3.gz")):
#         return "gff"
#     else:
#         sys.exit(
#             "ERROR: unrecognized extension of the annotation file.\n"
#             "Supported are gtf, gtf.gz, gff, gff.gz, gff2, gff2.gz, gff3 and gff3.gz"
#         )


# def transcripts_in_window(df, region, feature="transcript"):
#     """
#     Return the transcript names for which
#     either the end or the begin of an exon is within the window
#     """
#     return df.loc[
#         df["begin"].between(region.begin, region.end) | df["end"].between(
#             region.begin, region.end),
#         feature,
#     ].unique()


# def assign_colors_to_genes(transcripts):
#     genes = set([t.gene for t in transcripts])
#     colordict = {g: c for g, c in zip(genes, plcolors * 100)}
#     for t in transcripts:
#         t.color = colordict[t.gene]


# def parse_annotation(gtff, region, simplify=False):
#     type = annot_file_sniffer(gtff)
#     logging.info(f"Parsing {type} file...")
#     """
#     Parse the gtff and select the relevant region as determined by the window
#     return as Transcript objects
#     """
#     if Path(gtff + ".tbi").is_file():
#         logging.info(f"Reading {gtff} using a tabix stream.")
#         rg = f"{region.chromosome}:{region.begin}-{region.end}"
#         try:
#             tabix_stream = subprocess.Popen(['tabix', gtff, rg],
#                                             stdout=subprocess.PIPE,
#                                             stderr=subprocess.PIPE)
#         except FileNotFoundError as e:
#             logging.error("Error when opening a tabix stream.")
#             logging.error(e, exc_info=True)
#             sys.stderr.write(
#                 "\n\nERROR when opening a tabix stream.\n")
#             sys.stderr.write("Is tabix installed and on the PATH?.")
#             raise
#     annotationfile = pd.read_csv(tabix_stream.stdout, sep='\t', header=None, names=[
#                                  "chromosome", "source", "feature_type", "begin", "end", "score", "strand", "phase", "attributes"])
#     annotationfile["attributes"] = annotationfile.attributes.str.split(";")
#     for x in annotationfile.attributes:
#         for i in x:
#             if i.startswith("gene_name"):
#                 annotationfile['col1'] = i
#             if i.startswith("transcript_id"):
#                 annotationfile['col2'] = i
#             # if i.startswith("locus_tag"):
#             #     annotationfile["col3"] = i
#     annotationfile["gene"] = annotationfile.col1.str.split('=').str[-1]
#     annotationfile["transcript"] = annotationfile.col2.str.split('=').str[-1]
#     annotationfile.drop(columns=["source", "feature_type", "score",
#                         "phase", "attributes", "col1", "col2"], inplace=True)
#     if simplify:
#         annotationfile.drop_duplicates(
#             subset=["chromosome", "begin", "end", "gene"], inplace=True)
#         res = []
#         sys.stderr.write("0")
#         for g in transcripts_in_window(annotationfile, region, feature="gene"):
#             sys.stderr.write("1")
#             gtable = annotationfile.loc[annotationfile["gene"] == g]
#             sys.stderr.write("2")
#             if len(gtable):
#                 sys.stderr.write("3")
#                 res.append(
#                     Transcript(
#                         transcript=gtable["gene"].tolist()[0],
#                         gene=gtable["gene"].tolist()[0],
#                         exon_tuples=gtable.loc[:, ["begin", "end"]]
#                         .sort_values("begin")
#                         .itertuples(index=False, name=None),
#                         strand=gtable["strand"].tolist()[0],
#                     )
#                 )
#         sys.stderr.write(f"Found {len(res)} gene(s) in the region.\n")
#         logging.info(f"Found {len(res)} gene(s) in the region.\n")
#     else:
#         res = []
#         for t in transcripts_in_window(annotationfile, region, feature="transcript"):
#             tr = annotationfile.loc[annotationfile["transcript"] == t]
#             res.append(
#                 Transcript(
#                     transcript=t,
#                     gene=tr["gene"].tolist()[0],
#                     exon_tuples=tr.loc[:, ["begin", "end"]]
#                     .sort_values("begin")
#                     .itertuples(index=False, name=None),
#                     strand=tr["strand"].tolist()[0],
#                 )
#             )
#         sys.stderr.write(f"Found {len(res)} transcript(s) in the region.\n")
#         logging.info(f"Found {len(res)} transcript(s) in the region.\n")
#     assign_colors_to_genes(res)
#     return res


# def make_per_gene_annot_line_trace(transcript, region, y_pos):
#     """Generate a line trace for the gene
#     Trace can get limited by the window sizes
#     """
#     return go.Scatter(x=[max(transcript.begin, region.begin),
#                          min(transcript.end, region.end)],
#                       y=[y_pos, y_pos],
#                       mode='lines',
#                       line=dict(width=2, color=transcript.color),
#                       name=transcript.transcript,
#                       text=transcript.gene,
#                       hoverinfo='text',
#                       showlegend=False)


# def make_per_exon_arrow_trace(transcript, begin, end, y_pos):
#     """Generate a line+marker trace for the exon
#     The shape is an arrow, as defined by the strand in transcript.marker
#     """
#     return go.Scatter(x=[begin, end],
#                       y=[y_pos, y_pos],
#                       mode='lines+markers',
#                       line=dict(width=8, color=transcript.color),
#                       name=transcript.transcript,
#                       text=transcript.gene,
#                       hoverinfo='text',
#                       showlegend=False,
#                       marker=dict(symbol=transcript.marker,
#                                   size=8))


# joined, heatmap = args_heatmap(
#     args.files, args.window, args.gff, simplify=True)


# ######heatmap-DENDROGRAM


# def get_args():
#     parser = ArgumentParser(
#         description="Create heatmap of methylation frequency")
#     action = parser.add_mutually_exclusive_group(required=True)
#     action.add_argument("--files",
#                         nargs="+", help="Haplotype specific frequency files.")
#     action.add_argument(
#         "--joinedtable", help="Table with combined information")
#     parser.add_argument("-w", "--window",
#                         help="window (region) to which the visualisation has to be restricted",
#                         required=True if "--example" not in sys.argv else False)  # chr17:44345246-44353106
#     parser.add_argument("--names", nargs="+", default=[])
#     parser.add_argument("--outtable",
#                         help="File to write results to.")
#     parser.add_argument("--outfig",
#                         help="File to write results to.")

#     args = parser.parse_args()
#     if len(args.names) == 0:
#         for b in args.files:
#             c = b.rsplit('.', 2)[0]
#             args.names.append(os.path.basename(c))
#     if len(args.files) != len(args.names):
#         d = len(args.files)
#         e = len(args.names)
#         sys.exit(
#             "input and names should have same length, length files = %s and length names = %s" % (d, e))
#     return args


# args = get_args()


# class Region(object):
#     def __init__(self, region, fasta=None):
#         if ':' in region:
#             try:
#                 self.chromosome, interval = region.replace(',', '').split(':')
#                 try:
#                     # see if just integer chromosomes are used
#                     self.chromosome = int(self.chromosome)
#                 except ValueError:
#                     pass
#                 self.begin, self.end = [int(i) for i in interval.split('-')]
#                 self.start = self.begin
#             except ValueError:
#                 sys.exit("\n\nERROR: Window (-w/--window) inproperly formatted, "
#                          "examples of accepted formats are:\n"
#                          "'chr5:150200605-150423790' or 'ENST00000647408'\n\n")
#             self.size = self.end - self.begin
#             if not self.size > 0:
#                 sys.exit(
#                     "\n\nERROR: Window (-w/--window) inproperly formatted, "
#                     "begin of the interval has to be smaller than end\n\n"
#                 )
#             self.string = f"{self.chromosome}_{self.begin}_{self.end}"


# def args_heatmap(files, window):
#     dfs = []
#     region = Region(window)
#     for file, name in zip(files, args.names):
#         if window:
#             if Path(file + ".tbi").is_file():
#                 logging.info(f"Reading {file} using a tabix stream.")
#                 rg = f"{region.chromosome}:{region.begin}-{region.end}"
#                 try:
#                     tabix_stream = subprocess.Popen(['tabix', file, rg],
#                                                     stdout=subprocess.PIPE,
#                                                     stderr=subprocess.PIPE)
#                 except FileNotFoundError as e:
#                     logging.error("Error when opening a tabix stream.")
#                     logging.error(e, exc_info=True)
#                     sys.stderr.write(
#                         "\n\nERROR when opening a tabix stream.\n")
#                     sys.stderr.write("Is tabix installed and on the PATH?.")
#                     raise
#                 header = gzip.open(file, 'rt').readline().rstrip().split('\t')
#                 df = pd.read_csv(tabix_stream.stdout,
#                                  sep='\t', header=None, names=header)
#         else:
#             df = pd.read_csv(file, sep="\t")

#         df = pr.PyRanges(df.drop(['group_sequence', 'called_sites_methylated',
#                                   'num_motifs_in_group', "called_sites"], axis=1).rename(
#             columns={"start": "Start", "chromosome": "Chromosome", "end": "End"}))
#         logging.info("Read the file in a dataframe.")

#         table = df.df

#         table["MeanPosition"] = (table["Start"] + table["End"])/2

#         table = table.set_index("MeanPosition").drop(
#             ["Chromosome", "Start", "End"], axis=1)

#         dfs.append(table.rename(columns={'methylated_frequency': name}))
#     joined = dfs[0].join(dfs[1:], how='outer')

#     joined = joined.replace('na', np.NaN).interpolate(axis=1)
#     joined.to_csv(args.outtable, sep="\t", na_rep=np.NaN, header=True)
#     columns = list(joined.columns.values)
#     rows = list(joined.index)

#     heatmap = dash_bio.Clustergram(
#         data=joined.loc[rows].values, column_labels=columns, row_labels=rows, line_width=2)

#     html = heatmap.to_html()

#     with open(args.outfig, 'w') as f:
#         f.write(html)
#     return joined, heatmap


# joined, heatmap = args_heatmap(args.files, args.window)

# #####heatmap-ordered samples (case for gfpt2)


# def get_args():
#     parser = ArgumentParser(
#         description="Create heatmap of methylation frequency")
#     parser.add_argument(
#         "--table", help="Table with combined information")
#     parser.add_argument("--outfig",
#                         help="File to write results to.")
#     args = parser.parse_args()
#     return args


# args = get_args()


# def args_heatmap(table):
#     joined = pd.read_table(table).set_index("MeanPosition")
#     heatmap = px.imshow(joined, text_auto=True)
#     heatmap.update_xaxes(tickangle=45, tickfont=dict(size=8))
#     heatmap.update_yaxes(tickfont=dict(size=10))
#     heatmap.update_layout(
#         title="MEI1 region chr22:41,699,265-41,699,802: methylation frequency per position for each haplotype")
#     html = heatmap.to_html()

#     with open(args.outfig, 'w') as f:
#         f.write(html)
#     return heatmap


# heatmap = args_heatmap(args.table)
