import os
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser(description="Create heatmap of methylation frequencies.")
    parser.add_argument("--input", help="input file")
    parser.add_argument("--outputdir", help="output directory")
    return parser.parse_args()


args = get_args()

# Try reading the input data with error handling for empty files
try:
    df = pd.read_table(args.input, sep="\t")
except pd.errors.EmptyDataError:
    print(
        f"The input file {args.input} is empty. Creating empty output files for each chunk."
    )
    # Initialize an empty dataframe with column names if the file is empty
    df = pd.DataFrame(columns=["chrom", "start", "end"])

# Process data (only if not empty)
if not df.empty:
    df[["chrom", "start", "end"]] = df["position"].str.split(" ", expand=True)
    df.drop(["position", "start"], axis=1, inplace=True)
    df["end"] = df["end"].astype(
        int
    )  # modkit files are zero based, so to make it 1-based i take the end position

basename = os.path.basename(args.input)
basename_noextension = basename.replace(".tsv", "")

base_pair_chunk_size = 25000000
max_chunk_end = 300000000
current_chunk_start = 0
current_chunk_end = base_pair_chunk_size
chunk_index = 1

while current_chunk_start < max_chunk_end:
    # Ensure the current chunk end does not exceed the maximum allowed end
    current_chunk_end = min(current_chunk_end, max_chunk_end)

    # If df is empty, create an empty chunk like chunk_df.empty scenario
    if df.empty:
        chunk_df = pd.DataFrame([[np.nan] * len(df.columns)], columns=df.columns)
    else:
        # Find rows within the current chunk's base pair range
        chunk_df = df[
            (df["end"] >= current_chunk_start) & (df["end"] < current_chunk_end)
        ]
        if chunk_df.empty:
            chunk_df = pd.DataFrame([[np.nan] * len(df.columns)], columns=df.columns)

    # Define the filename
    filename = os.path.join(
        args.outputdir,
        f"{basename_noextension}_{current_chunk_start + 1}_{current_chunk_end}.tsv",
    )

    chunk_df.to_csv(filename, sep="\t", header=True, index=False)

    # Move to the next chunk
    current_chunk_start = current_chunk_end
    current_chunk_end += base_pair_chunk_size
    chunk_index += 1
