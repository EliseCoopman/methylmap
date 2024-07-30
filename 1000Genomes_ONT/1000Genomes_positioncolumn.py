import os
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser

def get_args():
    parser = ArgumentParser(description="Create heatmap of methylation frequencies.")
    parser.add_argument("--input", help="input file")
    return parser.parse_args()


args = get_args()


def add_position_column(file):
    df = pd.read_table(file, sep="\t")
    if df["chrom"].isna().all() and df["end"].isna().all():
        df["position"] = np.nan
    else:
        df["position"] = df["chrom"] + ":" + df["end"].astype(str)
    df.drop(["chrom", "end"], axis=1, inplace=True)

    df.to_csv(sys.stdout, sep="\t", index=False, header=True, na_rep=np.nan)


add_position_column(args.input)
