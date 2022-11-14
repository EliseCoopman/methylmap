import methylmap.plots as plots
import methylmap.annotation as annotation
from methylmap.import_data import read_mods
from methylmap.import_data import Region

import os
import sys
import itertools
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
        fasta=args.fasta
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
    parser.add_argument("--fasta", help="Fasta reference file")
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
    fasta=False
):
    if window:
        window = Region(window, expand)

    num_col = 2 if gff else 1  # number of subplots (columns) needed
    subplots = plots.create_subplots(num_col)

    # frequencies table with all meth frequencies of all samples
    meth_data, window = read_mods(
        files, table, names, window, groups, gff, outtable, fasta
    ) 
    fig = plots.plot_methylation(subplots, meth_data, num_col)

    if gff:
        annotation_traces = annotation.gff_annotation(gff, window, simplify)
        for annot_trace in annotation_traces:
            fig.append_trace(trace=annot_trace, row=1, col=1)
        fig.update_xaxes(title_text="", showticklabels=False, zeroline=False, row=1, col=1)
        fig.update_yaxes(title_text="", showticklabels=True, zeroline=False, row=1, col=1)
    with open(outfig, "w") as f:
        f.write(fig.to_html())


if __name__ == "__main__":
    main()
