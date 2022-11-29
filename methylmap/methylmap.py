import methylmap.plots as plots
import methylmap.annotation as annotation
from methylmap.import_data import read_mods
from methylmap.import_data import Region
from methylmap.version import __version__

import os
import sys
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
        mod=args.mod
    )


def get_args():
    parser = ArgumentParser(description="Plotting tool for population scale nucleotide modifications.")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument(
        "-f",
        "--files",
        nargs="+",
        help="Nanopolish calculate_methylation_frequency.py output or BAM/CRAM files.",
    )
    action.add_argument("-t", "--table", help="Methfrequencytable or overviewtable input.")
    parser.add_argument(
        "-w",
        "--window",
        help="Region to visualise. Format: chr:start-end (Example: chr20:58839718-58911192)",
    )
    parser.add_argument("-n", "--names", nargs="*", default=[],help="List with sample names.")
    parser.add_argument("--gff", "--gtf", help="Add annotation track based on GTF/GFF file.")
    parser.add_argument(
        "--expand",
        help="Number of base pairs to expand the window with in both directions.",
        type=int,
        default=0,
    )
    parser.add_argument("--outtable", help="File to write the frequencies table to.")
    parser.add_argument("--outfig", help="File to write output heatmap (in HTML format) to.")
    parser.add_argument("--groups", nargs="*", help="List of experimental group for each sample.")
    parser.add_argument("-s", "--simplify", action="store_true", help="Simplify annotation track to show genes rather than transcripts.")  # default: False
    parser.add_argument("--fasta", help="Fasta reference file, required when input is BAM/CRAM files or overviewtable with BAM/CRAM files.")
    parser.add_argument(
        "--mod",
        help="Modified base of interest when BAM/CRAM files as input. Options are: 5mC, 5hmC, 5fC, 5caC, 5hmU, 5fU, 5caU, 6mA, 5oxoG, Xao. Default = 5mC",
        default="5mC",
        choices=["5mC", "5hmC", "5fC", "5caC", "5hmU", "5fU", "5caU", "6mA", "5oxoG", "Xao"],
    )
    parser.add_argument(
        "-v",
        "--version",
        help="Print version and exit.",
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
    outfig="browser.html",
    outtable=False,
    names=False,
    window=False,
    expand=False,
    gff=False,
    groups=False,
    simplify=False,
    fasta=False,
    mod=False
):
    if window:
        window = Region(window, expand)

    num_col = 2 if gff else 1  # number of subplots (columns) needed
    subplots = plots.create_subplots(num_col)

    # frequencies table with all meth frequencies of all samples
    meth_data, window = read_mods(files, table, names, window, groups, gff, outtable, fasta, mod)
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
