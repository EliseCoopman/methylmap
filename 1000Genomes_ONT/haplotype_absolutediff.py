import os
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser(description="Calculate absolute difference of methylation status between two haplotypes.")
    parser.add_argument("--input", help="input file")
    parser.add_argument("--output", help="output file")
    return parser.parse_args()


args = get_args()

df = pd.read_table(args.input, sep="\t")
df = df.rename(columns={"end": "position"})
df.set_index(["chrom","position"], inplace=True)


# Extract column names
column_names = df.columns.tolist()

# Remove suffixes '_H1' or '_H2'
cleaned_column_names = [name.rsplit("_H", 1)[0] for name in column_names]

samples = list(dict.fromkeys(cleaned_column_names))

for sample in samples:
    h1_col = f'{sample}_H1'
    h2_col = f'{sample}_H2'
    result_col = f'{sample}'
    df[result_col] = abs(df[h1_col] - df[h2_col])
    df.drop(columns=[h1_col, h2_col], inplace=True)

df.to_csv(sys.std, sep="\t", na_rep=np.NaN, index=True)