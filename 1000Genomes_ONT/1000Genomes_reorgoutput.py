import os
import sys
import pandas as pd
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser(description="Create heatmap of methylation frequencies.")
    parser.add_argument("--input", help="input file")
    return parser.parse_args()


args = get_args()


def process_bed_files(input):
    columns_of_interest = [0, 1, 2, 10]  # Chromosome, start, end, methylation frequency

    bedfile = os.path.join(input)

    if not os.path.exists(bedfile):
        print(f"{bedfile} does not exist.")
        return

    df1 = pd.read_csv(bedfile, sep="\t", header=None)
    df1_selected = df1[columns_of_interest]
    df1_selected.to_csv(sys.stdout, sep="\t", header=False, index=False)

process_bed_files(args.input)
