import os
import sys
import numpy as np
import pandas as pd
import plotly.express as px
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser(description="Create plots.")
    parser.add_argument("--input", help="input file")
    parser.add_argument("--output_fig1", help="output file")
    parser.add_argument("--output_fig2", help="output file")
    return parser.parse_args()


args = get_args()

df = pd.read_table(args.input, sep="\t")
df.set_index(["chrom", "position"], inplace=True)

df_cleaned = df.dropna(how="all")
df_cleaned = df_cleaned[df_cleaned.index.get_level_values("chrom").isin(["chr1", "chr2", "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"])]


# Round the values in the DataFrame to 0 decimal places
df_rounded = df_cleaned.round(0)

# Flatten the rounded DataFrame into a Series
values = df_rounded.values.flatten()

# Count the frequency of each rounded value
value_counts = pd.Series(values).value_counts().sort_index()

# Create a DataFrame for Plotly
plot_data = pd.DataFrame({"Value": value_counts.index, "Count": value_counts.values})

# Create a bar plot using Plotly
fig = px.bar(
    plot_data,
    x="Value",
    y="Count",
    log_y=True,
    labels={"Value": "Values", "Count": "Count"},
    title="Occurrence of Rounded Values in DataFrame",
)

# Save the plot to a file
html = fig.to_html()
with open(args.output_fig1, "w") as f:
    f.write(html)


df_cleaned["percentage_90_or_above"] = df_cleaned.apply(
    lambda row: (
        ((row >= 90).sum() / (row.notna().sum())) * 100 if row.notna().sum() > 0 else 0
    ),
    axis=1,
)
df_cleaned["rounded_percentage"] = df_cleaned["percentage_90_or_above"].round(0)
counts = df_cleaned["rounded_percentage"].value_counts().sort_index().reset_index()
counts.columns = ["rounded_percentage", "count"]

fig2 = px.bar(
    counts,
    x="rounded_percentage",
    y="count",
    log_y=True,
    labels={
        "rounded_percentage": "Percentage of Values >= 90 (Rounded)",
        "count": "Count",
    },
    title="Distribution of Rounded Percentages of Values >= 90",
)

html = fig2.to_html()
with open(args.output_fig2, "w") as f:
    f.write(html)
