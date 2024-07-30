import os
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser

def get_args():
    parser = ArgumentParser(description="Create heatmap of methylation frequencies.")
    parser.add_argument("--input", nargs="+", help="input file")
    return parser.parse_args()


args = get_args()

dfs = []
for file in args.input:
    df = pd.read_table(
        file,
        sep="\t",
    )
    dfs.append(df.set_index("position"))

joined = dfs[0].join(dfs[1:], how="outer")
joined.reset_index(inplace=True)


if joined["position"].isna().sum() == 0:
    joined[["chr", "end"]] = df["position"].str.split(":", expand=True)
    joined["end"] = pd.to_numeric(joined["end"])
    joined.sort_values(by="end", inplace=True)
    joined.drop(columns=["chr", "end"], inplace=True)
else:
    non_na_rows = joined[joined["position"].notna()]
    na_rows = joined[joined["position"].isna()]

    if not non_na_rows.empty:
        non_na_rows[["chr", "end"]] = non_na_rows["position"].str.split(
            ":", expand=True
        )
        non_na_rows["end"] = pd.to_numeric(non_na_rows["end"])
        non_na_rows.sort_values(by="end", inplace=True)
        non_na_rows.drop(columns=["chr", "end"], inplace=True)
        joined = pd.concat([non_na_rows, na_rows], ignore_index=True)
    else:
        joined = na_rows

joined.set_index("position", inplace=True)
joined.to_csv(sys.stdout, sep="\t", na_rep=np.NaN, index=True)
