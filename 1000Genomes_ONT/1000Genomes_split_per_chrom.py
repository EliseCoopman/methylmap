import os
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser

def get_args():
    parser = ArgumentParser(description="Create heatmap of methylation frequencies.")
    parser.add_argument("--input", help="input file")
    parser.add_argument("--outputdir", help="output file")
    return parser.parse_args()


def process_file(file):
    df = pd.read_table(file, header=None, sep="\t")

    df["position"] = (
        df[0].astype(str) + " " + df[1].astype(str) + " " + df[2].astype(str)
    )
    basename = os.path.basename(file)
    basename_noextension = basename.replace(".bed.gz", "")
    new_column_names = ["chrom", "start", "end", basename_noextension]
    df.rename(
        columns={
            df.columns[i]: new_column_names[i] for i in range(len(new_column_names))
        },
        inplace=True,
    )
    return df.set_index("position"), basename_noextension


def split_by_chromosome(df):
    return {
        chrom: df_group.drop(["chrom", "start", "end"], axis=1)
        for chrom, df_group in df.groupby("chrom")
    }


args = get_args()


df, basename_noextension = process_file(args.input)
chromosome_data = split_by_chromosome(df)

for chrom, data in chromosome_data.items():
    output_file = os.path.join(args.outputdir, f"{basename_noextension}_{chrom}.tsv")
    data.to_csv(output_file, sep="\t", header=True, index=True)
