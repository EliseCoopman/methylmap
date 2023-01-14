import methylmap.plots as plots
import methylmap.annotation as annotation
import methylmap.dendro as dendrogram
from methylmap.import_data import read_mods
from methylmap.import_data import Region
from methylmap.version import __version__

import os
import sys
import numpy as np
from argparse import ArgumentParser


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
        fasta=args.fasta,
        mod=args.mod,
        dendro=args.dendro
    )


def get_args():
    parser = ArgumentParser(description="Plotting tool for population scale nucleotide modifications.")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument(
        "-f",
        "--files",
        nargs="+",
        help="list with nanopolish (processed with calculate_methylation_frequency.py) files or BAM/CRAM files",
    )
    action.add_argument("-t", "--table", help="methfreqtable or overviewtable input")
    parser.add_argument(
        "-w",
        "--window",
        help="region to visualise, format: chr:start-end (example: chr20:58839718-58911192)",
    )
    parser.add_argument("-n", "--names", nargs="*", default=[],help="list with sample names")
    parser.add_argument("--gff", "--gtf", help="add annotation track based on GTF/GFF file")
    parser.add_argument(
        "--expand",
        help="number of base pairs to expand the window with in both directions",
        type=int,
        default=0,
    )
    parser.add_argument("--outtable", help="file to write the frequencies table to in tsv format")
    parser.add_argument("--outfig", help="file to write output heatmap to, default: methylmap_{chr}_{start}_{end}.html (missing paths will be created)")
    parser.add_argument("--groups", nargs="*", help="list of experimental group for each sample")
    parser.add_argument("-s", "--simplify", action="store_true", help="simplify annotation track to show genes rather than transcripts")  # default: False
    parser.add_argument("--fasta", help="fasta reference file, required when input is BAM/CRAM files or overviewtable with BAM/CRAM files")
    parser.add_argument(
        "--mod",
        help="modified base of interest when BAM/CRAM files as input. Options are: 5mC, 5hmC, 5fC, 5caC, 5hmU, 5fU, 5caU, 6mA, 5oxoG, Xao, default = 5mC",
        default="5mC",
        choices=["5mC", "5hmC", "5fC", "5caC", "5hmU", "5fU", "5caU", "6mA", "5oxoG", "Xao"],
    )
    parser.add_argument("--dendro", action="store_true", help="perform hierarchical clustering on the samples/haplotypes and visualize with dendrogram on sorted heatmap as output")
    parser.add_argument(
        "-v",
        "--version",
        help="print version and exit",
        action="version",
        version=f"methylmap {__version__}",
    )
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
    outfig=None,
    outtable=False,
    names=False,
    window=False,
    expand=False,
    gff=False,
    groups=False,
    simplify=False,
    fasta=False,
    mod=False,
    dendro=False
):
    if window:
        window = Region(window, expand)

    num_col = 2 if gff else 1  # number of subplots (columns) needed
    num_row = 2 if dendro else 1
    subplots = plots.create_subplots(num_col, num_row)

    # frequencies table with all meth frequencies of all samples
    meth_data, window = read_mods(files, table, names, window, groups, gff, fasta, mod, dendro)
    if dendro:
        meth_data, den, list_sorted_samples = dendrogram.make_dendro(meth_data, window)
    meth_data.to_csv(outtable, sep="\t", na_rep=np.NaN, header=True)
    fig = plots.plot_methylation(subplots, meth_data, num_col, num_row)

    if gff:
        annotation_traces = annotation.gff_annotation(gff, window, simplify)
        for annot_trace in annotation_traces:
            fig.append_trace(trace=annot_trace, row=num_row, col=1)
        fig.update_xaxes(title_text="", showticklabels=False, zeroline=False, row=num_row, col=1)
        fig.update_yaxes(title_text="", showticklabels=True, zeroline=False, row=num_row, col=1)
    
    if dendro:
        for trace in den.select_traces():
            fig.add_trace(trace, row=1, col=num_col)
        fig.update_xaxes(title_text="", showticklabels=False, zeroline=False, showgrid=False, row=1, col=num_col)
        fig.update_yaxes(title_text="", showticklabels=False, zeroline=False, showgrid=False, row=1, col=num_col)
        fig.update_layout(showlegend=False)
        fig['data'][0]['x'] = den.layout.xaxis.tickvals

    if dendro and gff:
        fig['layout']['xaxis4']['tickvals'] = den.layout.xaxis.tickvals
        fig['layout']['xaxis4']['ticktext'] = list_sorted_samples
    if dendro and not gff:
        fig['layout']['xaxis2']['tickvals'] = den.layout.xaxis.tickvals
        fig['layout']['xaxis2']['ticktext'] = list_sorted_samples

    plots.create_output_methylmap(fig,outfig, window)

if __name__ == "__main__":
    main()
