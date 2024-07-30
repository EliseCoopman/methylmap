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
    df_cleaned = df.dropna(how="all")
    if not df_cleaned.empty:
        dfs.append(df_cleaned.set_index("position"))


result = pd.concat(dfs)
result.reset_index(inplace=True)
result[["chrom", "end"]] = result["position"].str.split(":", expand=True)
result["end"] = pd.to_numeric(result["end"])
result.sort_values(by="end", inplace=True)
result.drop(columns=["position"], inplace=True)
result = result[
    ["chrom", "end"] + [col for col in result.columns if col not in ["chrom", "end"]]
]

result.set_index(["chrom", "end"], inplace=True)
result.to_csv(sys.stdout, sep="\t", na_rep=np.NaN, index=True)
