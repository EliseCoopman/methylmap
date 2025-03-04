import methylmap.plots as plots
import methylmap.annotation as annot
import methylmap.dendro as dendrogram
from methylmap.import_data import read_mods
from methylmap.region import Region
from methylmap.process_1000Genomes import process_1000Genomes
from methylmap.version import __version__


import os
import re
import io
import sys
import json
import glob
import logging
import gzip
import pysam
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from argparse import ArgumentParser

import dash
import dash_bootstrap_components as dbc
from dash import Dash, dcc, html, Input, Output, State, ctx

app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler()],
    )
    args = get_args()
    annotation_dir = args.annotationdir
    file_paths = glob.glob(os.path.join(annotation_dir, "*_sorted.gff3.gz"))
    formatted_list = []
    for file_path in file_paths:
        file_name = os.path.basename(file_path).replace("_sorted.gff3.gz", "")
        entry = {
            "label": file_name.capitalize(),
            "value": file_name,
        }
        formatted_list.append(entry)

    db = args.db
    gff = args.gff
    genes_to_coords = make_gene_to_coords_dict(gff)

    button_style = {
        "height": "30px",
        "width": "100px",
        "background-color": "'#d3d3d3'",
        "color": "black",
        "align-items": "center",
        "vertical-align": "middle",
        "display": "inline-block",
        "justify-content": "center",
        "text-align": "center",
        "line-height": "20px",
        "border": "none",
        "font-size": "14px",
        "cursor": "pointer",
        "font-family": "Arial",
        "margin": "0px 2px",
    }
    tab_style = {
        "borderBottom": "1px solid #d6d6d6",
        "padding": "6px",
        "align-items": "center",
        "vertical-align": "middle",
        "display": "flex",
        "justify-content": "center",
        "align-items": "center",
    }

    tab_selected_style = {
        "borderTop": "1px solid #d6d6d6",
        "borderBottom": "1px solid #d6d6d6",
        "fontWeight": "bold",
        "padding": "6px",
        "align-items": "center",
        "vertical-align": "middle",
        "display": "flex",
        "justify-content": "center",
        "align-items": "center",
    }

    app.layout = html.Div(
        [
            html.Div(
                [
                    html.H1(
                        children="Methylmap: Visualization of modification frequencies for large cohort sizes",
                    ),
                ],
                style={
                    "textAlign": "center",
                },
            ),
            dcc.Tabs(
                [
                    dcc.Tab(
                        label="Overview",
                        style=tab_style,
                        selected_style=tab_selected_style,
                        children=[
                            dbc.Container(
                                [
                                    dbc.Row(
                                        [
                                            html.H1(
                                                "Example: GNAS region from ONT 1000Genomes",
                                                style={
                                                    "bottommargin": "0px",
                                                    "font-size": "25px",
                                                    "padding-left": "20px",
                                                },
                                            ),
                                        ],
                                    ),
                                    dbc.Row(
                                        [
                                            html.Img(
                                                src=dash.get_asset_url(
                                                    "1000Genomes_GNAS.png"
                                                ),
                                                style={"width": "100%"},
                                            ),
                                        ],
                                    ),
                                    dbc.Row(
                                        [
                                            html.Br(),
                                        ]
                                    ),
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                html.H1(
                                                    "About",
                                                    style={
                                                        "bottommargin": "0px",
                                                        "font-size": "25px",
                                                        "padding-left": "20px",
                                                    },
                                                ),
                                                width=2,
                                            ),
                                            dbc.Col(
                                                html.Div(
                                                    [
                                                        html.P(
                                                            [
                                                                "Welcome to the web application of methylmap, a tool for visualizing nucleotide modification frequencies of large cohort sizes. "
                                                                "In addition to the visualization of your own modification data, methylmap also provides the possibility to visualize the haplotype specific methylation frequencies of the ONT 1000Genomes dataset. "
                                                                "If this tool is useful for your research, please cite the following paper: ",
                                                                html.A(
                                                                    "Coopman et al.",
                                                                    href="https://www.biorxiv.org/content/10.1101/2022.11.28.518239v1",
                                                                    target="_blank",
                                                                ),
                                                                ", as well as the ONT 1000Genomes dataset: ",
                                                                html.A(
                                                                    "Gustafson et al.",
                                                                    href="https://www.medrxiv.org/content/10.1101/2024.03.05.24303792v1",
                                                                    target="_blank",
                                                                ),
                                                                ".",
                                                            ],
                                                            style={
                                                                "textAlign": "justify"
                                                            },
                                                        ),
                                                        html.P(
                                                            [
                                                                "Bam files of the 1000Genomes project were downloaded from ",
                                                                html.A(
                                                                    "hg38_bamfiles_1000GenomesONT",
                                                                    href="https://s3.amazonaws.com/1000g-ont/index.html?prefix=ALIGNMENT_AND_ASSEMBLY_DATA/FIRST_100/IN-HOUSE_MINIMAP2/HG38/",
                                                                    target="_blank",
                                                                ),
                                                                ", selected for only being base/modification called with dorado (5mCG and 5hmCG model), and further processed according to ",
                                                                html.A(
                                                                    "the pipeline",
                                                                    href="https://github.com/EliseCoopman/methylmap/blob/main/1000Genomes/1000Genomes_snakemake.smk",
                                                                ),
                                                                " available at the methylmap Github page. ",
                                                                "In total, 226 samples were processed, resulting in visualization of 452 haplotypes that can be consulted in the 'ONT 1000Genomes' tab.",
                                                            ],
                                                            style={
                                                                "textAlign": "justify"
                                                            },
                                                        ),
                                                        html.P(
                                                            [
                                                                "Your own data can be visualized in the 'Upload your own data' tab' by providing a modification frequency table as .tsv file. The maximum size of the file is 100MB. "
                                                                "Starting from bam or cram files with MM/ML tags, such a table for your genomic region of interest can be generated using the ",
                                                                html.A(
                                                                    "multiparsetable.py script",
                                                                    href="https://github.com/EliseCoopman/methylmap/blob/main/multiparsetable.py",
                                                                    target="_blank",
                                                                ),
                                                                " available at the methylmap Github page. This script also supports the generation of a modification frequency table from Nanopolish input files. "
                                                                "More information about the expected formating for the modification frequency table can be found at ",
                                                                html.A(
                                                                    "the methylmap github page",
                                                                    href="https://github.com/EliseCoopman/methylmap",
                                                                    target="_blank",
                                                                ),
                                                                ". "
                                                                "Any remarks, issues, suggestions or questions can be uploaded in the form of an ",
                                                                html.A(
                                                                    "issue on the github page.",
                                                                    href="https://github.com/EliseCoopman/methylmap/issues",
                                                                    target="_blank",
                                                                ),
                                                            ],
                                                            style={
                                                                "textAlign": "justify"
                                                            },
                                                        ),
                                                    ]
                                                ),
                                                width=8,
                                            ),
                                        ],
                                    ),
                                ],
                            ),
                        ],
                    ),
                    dcc.Tab(
                        label="ONT 1000Genomes",
                        style=tab_style,
                        selected_style=tab_selected_style,
                        children=[
                            html.H1(
                                "Haplotype-specific methylation frequencies ONT 1000Genomes",
                                style={
                                    "bottommargin": "0px",
                                    "font-size": "25px",
                                    "padding-left": "20px",
                                },
                            ),
                            dbc.Container(
                                children=[
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                html.Label("Genomic region/Gene name:"),
                                                width=3,
                                            ),
                                            dbc.Col(
                                                html.Div(
                                                    [
                                                        input_box_genomebrowser(
                                                            app, genes_to_coords, args
                                                        ),
                                                        html.Button(
                                                            id="confirm-button_1000Genomes",
                                                            n_clicks=0,
                                                            children="Confirm",
                                                            style={
                                                                "background-color": "'#A9A9A9'",
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "display": "flex",
                                                        "align-items": "center",
                                                    },
                                                )
                                            ),
                                            dbc.Col(
                                                html.Div(
                                                    id="error-message_1000Genomes_inputbox",
                                                    style={"color": "red"},
                                                ),
                                            ),
                                            dbc.Col(
                                                html.Div(
                                                    id="error-message_1000Genomes_dccstore",
                                                    style={"color": "red"},
                                                ),
                                            ),
                                        ],
                                    ),
                                    dbc.Row(style={"margin-top": "10px"}),
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                html.Label("Zoom in/out:"), width=3
                                            ),
                                            dbc.Col(
                                                html.Div(
                                                    [
                                                        html.Button(
                                                            id="button-o10_1000Genomes",
                                                            n_clicks=0,
                                                            children="Out 10x",
                                                            style=button_style,
                                                        ),
                                                        html.Button(
                                                            id="button-o3_1000Genomes",
                                                            n_clicks=0,
                                                            children="Out 3x",
                                                            style=button_style,
                                                        ),
                                                        html.Button(
                                                            id="button-i3_1000Genomes",
                                                            n_clicks=0,
                                                            children="In 3x",
                                                            style=button_style,
                                                        ),
                                                        html.Button(
                                                            id="button-i10_1000Genomes",
                                                            n_clicks=0,
                                                            children="In 10x",
                                                            style=button_style,
                                                        ),
                                                    ]
                                                )
                                            ),
                                        ],
                                    ),
                                    dbc.Row(
                                        [
                                            html.Br(),
                                        ]
                                    ),
                                    dbc.Row(
                                        [
                                            dbc.Col(html.Label("Annotation:"), width=3),
                                            dbc.Col(
                                                dcc.Dropdown(
                                                    id="annotation_1000Genomes",
                                                    options=[
                                                        {
                                                            "label": "On",
                                                            "value": "on",
                                                        },
                                                        {
                                                            "label": "Off",
                                                            "value": "off",
                                                        },
                                                    ],
                                                    value="on",
                                                    clearable=False,
                                                ),
                                                width=2,
                                            ),
                                        ],
                                        align="center",
                                    ),
                                    html.Div(
                                        id="annotation-type-container_1000Genomes",
                                        children=[
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        html.Label("Annotation type:"),
                                                        width=3,
                                                    ),
                                                    dbc.Col(
                                                        dcc.Dropdown(
                                                            id="annotation-type_1000Genomes",
                                                            options=[
                                                                {
                                                                    "label": "Gene",
                                                                    "value": "gene",
                                                                },
                                                                {
                                                                    "label": "Transcript",
                                                                    "value": "transcript",
                                                                },
                                                            ],
                                                            value="gene",
                                                            clearable=False,
                                                        ),
                                                        width=2,
                                                    ),
                                                ],
                                                align="center",
                                            ),
                                        ],
                                    ),
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                html.Label("Hierarchical clustering:"),
                                                width=3,
                                            ),
                                            dbc.Col(
                                                dcc.Dropdown(
                                                    id="hierarchical_clustering_1000Genomes",
                                                    options=[
                                                        {
                                                            "label": "On",
                                                            "value": "on",
                                                        },
                                                        {
                                                            "label": "Off",
                                                            "value": "off",
                                                        },
                                                    ],
                                                    value="off",
                                                    clearable=False,
                                                ),
                                                width=2,
                                            ),
                                            dbc.Col(
                                                html.Div(
                                                    id="hierarchical_clustering_1000Genomes_message",
                                                    children="Individuals with 40% or more missing data are removed. If there are still missing values, they are estimated using interpolation. Any remaining missing values after interpolation result in the removal of those individuals.",
                                                    style={"fontSize": "5px"},
                                                ),
                                            ),
                                        ],
                                        align="center",
                                    ),
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                html.Label("Color scale:"),
                                                width=3,
                                            ),
                                            dbc.Col(
                                                dcc.Dropdown(
                                                    id="color_scale_1000Genomes",
                                                    options=[
                                                        {
                                                            "label": "BlueRed",
                                                            "value": "bluered",
                                                        },
                                                        {
                                                            "label": "Greys",
                                                            "value": "greys",
                                                        },
                                                        {
                                                            "label": "Plasma",
                                                            "value": "plasma",
                                                        },
                                                    ],
                                                    value="plasma",
                                                    clearable=False,
                                                ),
                                                width=2,
                                            ),
                                        ],
                                    ),
                                ],
                            ),
                            html.Div(
                                id="error-message_1000Genomes", style={"color": "red"}
                            ),
                            Genome_browser(args, app, gff, genes_to_coords),
                            dcc_store_genome_browser(app, db, genes_to_coords, args),
                        ],
                    ),
                    dcc.Tab(
                        label="Upload your own data",
                        style=tab_style,
                        selected_style=tab_selected_style,
                        children=[
                            html.H1(
                                "Upload your own data",
                                style={
                                    "bottommargin": "0px",
                                    "font-size": "25px",
                                    "padding-left": "20px",
                                },
                            ),
                            dcc.Upload(
                                id="upload-data",
                                children=html.Div(
                                    [
                                        html.P(
                                            [
                                                "Drag and drop your .tsv/.tsv.gz modification frequency table, enter your genomic region of interest and click on the 'Confirm' button. Please note that the maximum allowed file size is 100MB."
                                            ],
                                            style={
                                                "textAlign": "center",
                                                "fontWeight": "bold",
                                                "textSize": "10px",
                                            },
                                        ),
                                        html.P(
                                            [
                                                "Starting from bam or cram files with MM/ML tags or Nanopolish files, such a table can be generated using the ",
                                                html.A(
                                                    "multiparsetable.py script",
                                                    href="https://github.com/EliseCoopman/methylmap/blob/main/multiparsetable.py",
                                                ),
                                                " available at the methylmap Github page. ",
                                                "This table should contain the modification frequencies for your genomic region of interest and is expected to have the following columns: chrom, position, samplename1, samplename2, ... .",
                                                " Not more than one unique chromosome is allowed in the table. For more information, please consult the ",
                                                html.A(
                                                    "methylmap Github page",
                                                    href="https://github.com/EliseCoopman/methylmap",
                                                ),
                                                ".",
                                            ],
                                            style={
                                                "textAlign": "justify",
                                                "textSize": "3px",
                                                "fontStyle": "italic",
                                                "margin": "8px",
                                            },
                                        ),
                                    ]
                                ),
                                style={
                                    "width": "100%",
                                    "height": "100px",
                                    "lineHeight": "60px",
                                    "borderWidth": "1px",
                                    "borderStyle": "dashed",
                                    "borderRadius": "5px",
                                    "line-height": "20px",
                                    "margin": "10px",
                                },
                                multiple=False,
                                max_size=100000000,
                            ),
                            dcc.Loading(
                                id="upload-loading",
                                type="default",
                                children=html.Div(
                                    id="loading-content",
                                ),
                            ),
                            dbc.Container(
                                children=[
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                html.Label("Genomic region:"), width=3
                                            ),
                                            dbc.Col(
                                                html.Div(
                                                    [
                                                        input_box(app, args),
                                                        html.Button(
                                                            id="confirm-button",
                                                            n_clicks=0,
                                                            children="Confirm",
                                                            style={
                                                                "background-color": "'#A9A9A9'",
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "display": "flex",
                                                        "align-items": "center",
                                                    },
                                                )
                                            ),
                                            dbc.Col(
                                                html.Div(
                                                    id="error-message-inputbox",
                                                    style={"color": "red"},
                                                ),
                                            ),
                                            dbc.Col(
                                                html.Div(
                                                    id="error-message-uploadowndata",
                                                    style={"color": "red"},
                                                ),
                                            ),
                                        ],
                                    ),
                                    dbc.Row(style={"margin-top": "10px"}),
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                html.Label("Zoom in/out:"), width=3
                                            ),
                                            dbc.Col(
                                                html.Div(
                                                    [
                                                        html.Button(
                                                            id="button-o10",
                                                            n_clicks=0,
                                                            children="Out 10x",
                                                            style=button_style,
                                                        ),
                                                        html.Button(
                                                            id="button-o3",
                                                            n_clicks=0,
                                                            children="Out 3x",
                                                            style=button_style,
                                                        ),
                                                        html.Button(
                                                            id="button-i3",
                                                            n_clicks=0,
                                                            children="In 3x",
                                                            style=button_style,
                                                        ),
                                                        html.Button(
                                                            id="button-i10",
                                                            n_clicks=0,
                                                            children="In 10x",
                                                            style=button_style,
                                                        ),
                                                    ]
                                                )
                                            ),
                                        ],
                                    ),
                                    dbc.Row(
                                        [
                                            html.Br(),
                                        ]
                                    ),
                                    dbc.Row(
                                        [
                                            dbc.Col(html.Label("Annotation:"), width=3),
                                            dbc.Col(
                                                dcc.Dropdown(
                                                    id="annotation",
                                                    options=[
                                                        {
                                                            "label": "On",
                                                            "value": "on",
                                                        },
                                                        {
                                                            "label": "Off",
                                                            "value": "off",
                                                        },
                                                    ],
                                                    value="off",
                                                    clearable=False,
                                                ),
                                                width=5,
                                            ),
                                        ],
                                        align="center",
                                    ),
                                    html.Div(
                                        id="annotation-file-container",
                                        children=[
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        html.Label("Annotation file:"),
                                                        width=3,
                                                    ),
                                                    dbc.Col(
                                                        dcc.Dropdown(
                                                            id="annotation-file",
                                                            options=formatted_list,
                                                            clearable=False,
                                                        ),
                                                        width=5,
                                                    ),
                                                ],
                                                align="center",
                                            ),
                                        ],
                                    ),
                                    html.Div(
                                        id="annotation-type-container",
                                        children=[
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        html.Label("Annotation type:"),
                                                        width=3,
                                                    ),
                                                    dbc.Col(
                                                        dcc.Dropdown(
                                                            id="annotation-type",
                                                            options=[
                                                                {
                                                                    "label": "Gene",
                                                                    "value": "gene",
                                                                },
                                                                {
                                                                    "label": "Transcript",
                                                                    "value": "transcript",
                                                                },
                                                            ],
                                                            value="gene",
                                                            clearable=False,
                                                        ),
                                                        width=5,
                                                    ),
                                                ],
                                                align="center",
                                            ),
                                        ],
                                    ),
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                html.Label("Hierarchical clustering:"),
                                                width=3,
                                            ),
                                            dbc.Col(
                                                dcc.Dropdown(
                                                    id="hierarchical_clustering",
                                                    options=[
                                                        {
                                                            "label": "On",
                                                            "value": "on",
                                                        },
                                                        {
                                                            "label": "Off",
                                                            "value": "off",
                                                        },
                                                    ],
                                                    value="off",
                                                    clearable=False,
                                                ),
                                                width=5,
                                            ),
                                            dbc.Col(
                                                html.Div(
                                                    id="hierarchical_clustering_message",
                                                    children="Individuals with 40% or more missing data are removed. If there are still missing values, they are estimated using interpolation. Any remaining missing values after interpolation result in the removal of those individuals.",
                                                    style={"fontSize": "5px"},
                                                ),
                                            ),
                                        ],
                                        align="center",
                                    ),
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                html.Label("Color scale:"),
                                                width=3,
                                            ),
                                            dbc.Col(
                                                dcc.Dropdown(
                                                    id="color_scale",
                                                    options=[
                                                        {
                                                            "label": "BlueRed",
                                                            "value": "bluered",
                                                        },
                                                        {
                                                            "label": "Greys",
                                                            "value": "greys",
                                                        },
                                                        {
                                                            "label": "Plasma",
                                                            "value": "plasma",
                                                        },
                                                    ],
                                                    value="plasma",
                                                    clearable=False,
                                                ),
                                                width=5,
                                            ),
                                        ],
                                    ),
                                ],
                            ),
                            html.Div(
                                id="error-message",
                                style={"color": "red"},
                            ),
                            meth_browser(app, args, gff, annotation_dir),
                            dcc_store(app, args),
                        ],
                    ),
                ]
            ),
        ],
    )

    @app.callback(
        Output("annotation-type-container", "style"), Input("annotation", "value")
    )
    def toggle_annotation_type_visibility(annotation):
        if "on" in annotation:
            return {"display": "block"}
        else:
            return {"display": "none"}

    @app.callback(
        Output("annotation-type-container_1000Genomes", "style"),
        Input("annotation_1000Genomes", "value"),
    )
    def toggle_annotation_type_visibility(annotation_1000Genomes):
        if "on" in annotation_1000Genomes:
            return {"display": "block"}
        else:
            return {"display": "none"}

    @app.callback(
        Output("annotation-file-container", "style"), Input("annotation", "value")
    )
    def toggle_annotation_file_visibility(annotation):
        if "on" in annotation:
            return {"display": "block"}
        else:
            return {"display": "none"}

    @app.callback(
        Output("hierarchical_clustering_1000Genomes_message", "style"),
        Input("hierarchical_clustering_1000Genomes", "value"),
    )
    def update_message_visibility(selected_value):
        if selected_value == "on":
            return {"display": "block"}
        return {"display": "none"}

    @app.callback(
        Output("hierarchical_clustering_message", "style"),
        Input("hierarchical_clustering", "value"),
    )
    def update_message_visibility(selected_value):
        if selected_value == "on":
            return {"display": "block"}
        return {"display": "none"}

    @app.callback(
        Output("upload-loading", "children"),
        Input("upload-data", "contents"),
        State("upload-data", "filename"),
    )
    def process_upload(contents, filename):
        print(filename)
        if contents is None:
            return [html.Div()]

        if not (filename.endswith('.tsv') or filename.endswith('.tsv.gz')):
            return html.Div(
            [
                html.H5(
                    f"Error: Only .tsv or .tsv.gz files are allowed. The file '{filename}' has an invalid extension.",
                    style={"fontSize": "12px", "margin": "8px", "color": "red"},
                ),
            ]
        )

        import base64
        content_type, content_string = contents.split(",")
        decoded_content = base64.b64decode(content_string)
        
        try:
            if filename.endswith(".gz"):
                with gzip.open(io.BytesIO(decoded_content), "rt", encoding="utf-8") as f:
                    header = f.readline().strip()
            else:
                with io.StringIO(decoded_content.decode("utf-8")) as f:
                    header = f.readline().strip()
                    print(header)
            if 'chrom' not in header or 'position' not in header:
                print("chrom or position not in header")
                return html.Div(
                [
                    html.H5(
                        f"Error: The file '{filename}' must have 'chrom' and 'position' columns in the first row.",
                        style={"fontSize": "12px", "margin": "8px", "color": "red"},
                    ),
                ]
            )

            return html.Div(
            [
                html.H5(
                    f"File '{filename}' uploaded successfully. Please enter a genomic region and click on the 'Confirm' button.",
                    style={
                        "fontSize": "12px",
                        "margin": "8px",
                    },
                ),
            ]
        )
        except Exception as e:
            return html.Div(
            [
                html.H5(
                    f"Error: Could not process the file '{filename}'. Please ensure it's a valid .tsv or .tsv.gz file.",
                    style={"fontSize": "12px", "margin": "8px", "color": "red"},
                ),
            ]
        )

    app.title = "methylmap"
    app.run(host=args.host, port=args.port, debug=args.debug)


def get_args():
    parser = ArgumentParser(
        description="Plotting tool for population scale nucleotide modifications."
    )
    action = parser.add_mutually_exclusive_group()
    action.add_argument(
        "-f",
        "--files",
        nargs="+",
        help="list with BAM/CRAM files or nanopolish (processed with calculate_methylation_frequency.py) files",
    )
    action.add_argument("-t", "--table", help="methfreqtable or overviewtable input")
    parser.add_argument(
        "-w",
        "--window",
        help="region to visualise, format: chr:start-end (example: chr20:58839718-58911192)",
    )
    parser.add_argument(
        "-n", "--names", nargs="*", default=[], help="list with sample names"
    )
    parser.add_argument(
        "--gff",
        help="add annotation track based on GFF3 file",
        required=True,
    )
    parser.add_argument("--output", help="TSV file to write the frequencies to."),
    parser.add_argument(
        "--groups", nargs="*", help="list of experimental group for each sample"
    )
    parser.add_argument(
        "-s",
        "--simplify",
        action="store_true",
        help="simplify annotation track to show genes rather than transcripts",
    )  # default: False
    parser.add_argument(
        "--fasta",
        help="fasta reference file, required when input is BAM/CRAM files or overviewtable with BAM/CRAM files",
    )
    parser.add_argument(
        "--mod",
        help="modified base of interest when BAM/CRAM files as input. Options are: m, h, default = m",
        default="m",
        choices=["m", "h"],  # methylation  # hydroxymethylation
    )
    parser.add_argument(
        "--hapl",
        action="store_true",
        help="display modification frequencies in input BAM/CRAM file for each haplotype (alternating haplotypes in methylmap)",
    )
    parser.add_argument(
        "--dendro",
        action="store_true",
        help="perform hierarchical clustering on the samples/haplotypes and visualize with dendrogram on sorted heatmap as output",
    )
    parser.add_argument(
        "--threads",
        help="number of threads to use when processing BAM/CRAM files",
        type=int,
        default=12,
    )
    parser.add_argument("--db", help="use 1000Genomes data", required=True)
    parser.add_argument(
        "--annotationdir",
        required=True,
        help="directory with annotation files",
    )
    parser.add_argument("--quiet", action="store_true", help="suppress modkit output")
    parser.add_argument(
        "--debug",
        help="Run the app in debug mode",
        action="store_true",
    )
    parser.add_argument(
        "--host", help="Host IP used to serve the application", default="127.0.0.1"
    )
    parser.add_argument(
        "--port", help="Port used to serve the application", default=8050, type=int
    )
    parser.add_argument(
        "-v",
        "--version",
        help="print version and exit",
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
    if args.files or args.table:
        if not args.window:
            sys.exit("ERROR: please provide a genomic region with --window")
    return args


def Genome_browser(args, app, gff, genes_to_coords):
    gnas = "chr20:58839718-58911192"

    @app.callback(
        [
            Output(component_id="plot_1000Genomes", component_property="children"),
            Output(
                component_id="error-message_1000Genomes", component_property="children"
            ),
        ],
        [
            Input(
                component_id="confirm-button_1000Genomes", component_property="n_clicks"
            ),
            Input(component_id="button-o3_1000Genomes", component_property="n_clicks"),
            Input(component_id="button-o10_1000Genomes", component_property="n_clicks"),
            Input(component_id="button-i3_1000Genomes", component_property="n_clicks"),
            Input(component_id="button-i10_1000Genomes", component_property="n_clicks"),
            Input(
                component_id="hierarchical_clustering_1000Genomes",
                component_property="value",
            ),
            Input(component_id="annotation_1000Genomes", component_property="value"),
            Input(
                component_id="annotation-type_1000Genomes", component_property="value"
            ),
            Input(component_id="color_scale_1000Genomes", component_property="value"),
            Input(
                component_id="intermediate-data_1000Genomes", component_property="data"
            ),
        ],
        [
            State(component_id="input-box_1000Genomes", component_property="children"),
        ],
    )
    def update_meth_browser(
        button_confirm_1000Genomes,
        button_o3_1000Genomes,
        button_o10_1000Genomes,
        button_i3_1000Genomes,
        button_i10_1000Genomes,
        hierarchical_clustering_1000Genomes,
        annotation_1000Genomes,
        annotation_type_1000Genomes,
        color_scale,
        mod_data_1000Genomes,
        inputbox_1000Genomes,
    ):
        window, dendro, annotation, simplify, num_row, num_col, subplots = (
            browser1000Genomes_information(
                args,
                button_confirm_1000Genomes,
                button_o3_1000Genomes,
                button_o10_1000Genomes,
                button_i3_1000Genomes,
                button_i10_1000Genomes,
                inputbox_1000Genomes,
                hierarchical_clustering_1000Genomes,
                annotation_1000Genomes,
                annotation_type_1000Genomes,
                gnas,
                genes_to_coords,
            )
        )
        if window is None:
            return None, None
        json_data_1000Genomes = json.loads(mod_data_1000Genomes)
        mod_data_1000Genomes = pd.DataFrame(
            json_data_1000Genomes["data"], columns=json_data_1000Genomes["columns"]
        )
        mod_data_1000Genomes.index = json_data_1000Genomes["index"]
        fig = process_fig(
            mod_data_1000Genomes,
            dendro,
            subplots,
            num_row,
            num_col,
            window,
            annotation,
            simplify,
            gff,
            color_scale,
        )
        return (
            html.Div(
                dcc.Graph(
                    figure=fig,
                    config={
                        "toImageButtonOptions": {
                            "format": "png",
                            "filename": "1000Genomesplot",
                            "height": 800,
                            "width": 1200,
                            "scale": 12,
                        }
                    },
                ),
                id="plot_1000Genomes",
            ),
            None,
        )

    return html.Div(id="plot_1000Genomes")


def dcc_store_genome_browser(app, db, genes_to_coords, args):
    gnas = "chr20:58839718-58911192"

    @app.callback(
        Output(component_id="intermediate-data_1000Genomes", component_property="data"),
        Output(
            component_id="error-message_1000Genomes_dccstore",
            component_property="children",
        ),
        [
            Input(
                component_id="confirm-button_1000Genomes", component_property="n_clicks"
            ),
            Input(component_id="button-o3_1000Genomes", component_property="n_clicks"),
            Input(component_id="button-o10_1000Genomes", component_property="n_clicks"),
            Input(component_id="button-i3_1000Genomes", component_property="n_clicks"),
            Input(component_id="button-i10_1000Genomes", component_property="n_clicks"),
            Input(component_id="input-box_1000Genomes", component_property="children"),
        ],
        [
            State(component_id="input-box_1000Genomes", component_property="children"),
        ],
    )
    def generate_data(
        button_confirm_1000Genomes,
        button_o3_1000Genomes,
        button_o10_1000Genomes,
        button_i3_1000Genomes,
        button_i10_1000Genomes,
        input_box_1000Genomes,
        input_box_1000Genomes2,
    ):
        window = window_input_1000Genomes(
            args, input_box_1000Genomes, genes_to_coords, gnas
        )
        if window == (
            "error",
            "Input not recognized. Please enter genomic region in a valid format. Example chr20:58,839,718-58,911,192 or chr20:58839718-58911192",
        ):
            return (
                "error",
                "Input not recognized. Please enter genomic region in a valid format. Example chr20:58,839,718-58,911,192 or chr20:58839718-58911192",
            )
        if window == ("The start position must be less than the end position."):
            return (None, "The start position must be lower than the end position.")
        if window == (
            "The region is too large. Please enter a region smaller than 1,000,000 bp."
        ):
            return (
                None,
                "The region is too large. Please enter a region smaller than 1,000,000 bp.",
            )
        if (
            window
            == "Invalid chromosome. Chromsome not recognized/not present in the data."
        ):
            return (
                None,
                "Invalid chromosome. Chromsome not recognized/not present in the data.",
            )
        if window == "No data found for the given interval.":
            return (None, "No data found for the given interval.")
        window_region = Region(window)
        mod_data_1000Genomes = process_1000Genomes(db, window_region)
        mod_data_1000Genomes = mod_data_1000Genomes.reset_index(
            level="chrom", drop=True
        )
        json_data_1000Genomes = mod_data_1000Genomes.to_json(orient="split")
        return json_data_1000Genomes, None

    return dcc.Store(id="intermediate-data_1000Genomes")


def input_box_genomebrowser(app, genes_to_coords, args):
    gnas_region = "chr20:58,839,718-58,911,192"

    @app.callback(
        Output(component_id="input-box_1000Genomes", component_property="children"),
        Output(
            component_id="error-message_1000Genomes_inputbox",
            component_property="children",
        ),
        [
            Input(
                component_id="confirm-button_1000Genomes", component_property="n_clicks"
            ),
            Input(component_id="button-o3_1000Genomes", component_property="n_clicks"),
            Input(component_id="button-o10_1000Genomes", component_property="n_clicks"),
            Input(component_id="button-i3_1000Genomes", component_property="n_clicks"),
            Input(component_id="button-i10_1000Genomes", component_property="n_clicks"),
        ],
        [State(component_id="input-box_1000Genomes", component_property="children")],
    )
    def update_value(
        button_confirm_1000Genomes,
        button_o3_1000Genomes,
        button_o10_1000Genomes,
        button_i3_1000Genomes,
        button_i10_1000Genomes,
        window,
    ):
        window = window["props"]["value"] if window else gnas_region
        if not validate_input_1000Genomes(window, args):
            window = window.upper()
            coords = genes_to_coords.get(window)
            if not coords:
                window = "error"
                return (
                    html.Div(
                        dcc.Input(type="text", value=window),
                        id="input-box_1000Genomes",
                        style={"height": "30px", "margin": "0px 2px"},
                    ),
                    None,
                )
            else:
                window = coords
        window = Region(window)
        if "button-o3_1000Genomes" == ctx.triggered_id:
            window = window * 3
        elif "button-o10_1000Genomes" == ctx.triggered_id:
            window = window * 10
        elif "button-i3_1000Genomes" == ctx.triggered_id:
            window = window / 3
        elif "button-i10_1000Genomes" == ctx.triggered_id:
            window = window / 10
        window = window.fmt
        return (
            html.Div(
                dcc.Input(type="text", value=window),
                id="input-box_1000Genomes",
                style={"height": "30px", "margin": "0px 2px"},
            ),
            None,
        )

    return html.Div(
        dcc.Input(type="text", value=gnas_region),
        id="input-box_1000Genomes",
        style={"height": "30px", "margin": "0px 2px"},
    )


def process_fig(
    mod_data,
    dendro,
    subplots,
    num_row,
    num_col,
    window,
    annotation,
    simplify,
    gff,
    color_scale,
    output=False,
):
    if dendro:
        logging.info("Performing hierarchical clustering")
        mod_data, den, list_sorted_samples = dendrogram.make_dendro(mod_data)
    if output is not False:
        mod_data.to_csv(output, sep="\t", na_rep=np.NaN, header=True)
    fig = plots.plot_methylation(subplots, mod_data, num_col, num_row, color_scale)
    if annotation:
        window = Region(window)
        annotation_traces = annot.gff_annotation(gff, window, simplify)
        for annot_trace in annotation_traces:
            fig.append_trace(trace=annot_trace, row=num_row, col=1)
        fig.update_xaxes(
            title_text="", showticklabels=False, zeroline=False, row=num_row, col=1
        )
        fig.update_yaxes(
            title_text="", showticklabels=True, zeroline=False, row=num_row, col=1
        )
    if dendro:
        for trace in den.select_traces():
            fig.add_trace(trace, row=1, col=num_col)
        fig.update_xaxes(
            title_text="",
            showticklabels=False,
            zeroline=False,
            showgrid=False,
            row=1,
            col=num_col,
        )
        fig.update_yaxes(
            title_text="",
            showticklabels=False,
            zeroline=False,
            showgrid=False,
            row=1,
            col=num_col,
        )
        fig.update_layout(showlegend=False)
        fig["data"][0]["x"] = den.layout.xaxis.tickvals

        dendro_xaxis = "xaxis4" if annotation else "xaxis2"
        fig["layout"][dendro_xaxis]["tickvals"] = den.layout.xaxis.tickvals
        fig["layout"][dendro_xaxis]["ticktext"] = list_sorted_samples
    return fig


def browser1000Genomes_information(
    args,
    button_confirm,
    button_o3,
    button_o10,
    button_i3,
    button_i10,
    input_box,
    hierarchical_clustering,
    annotation,
    annotation_type,
    windowregion,
    genes_to_coords,
):
    if windowregion is None and input_box["props"]["value"] is None:
        return None, None, None, None, None, None, None
    window = window_input_1000Genomes(args, input_box, genes_to_coords, windowregion)
    if window == (
        "error",
        "Input not recognized. Please enter genomic region in a valid format. Example chr20:58,839,718-58,911,192 or chr20:58839718-58911192",
    ):
        return None, None, None, None, None, None, None
    if window == ("The start position must be less than the end position."):
        return None, None, None, None, None, None, None
    if window == (
        "The region is too large. Please enter a region smaller than 1,000,000 bp."
    ):
        return None, None, None, None, None, None, None
    if (
        window
        == "Invalid chromosome. Chromsome not recognized/not present in the data."
    ):
        return None, None, None, None, None, None, None
    if window == "No data found for the given interval.":
        return None, None, None, None, None, None, None
    dendro = "on" in hierarchical_clustering
    annotation = "on" in annotation
    if annotation:
        if annotation_type == "gene":
            simplify = True
        if annotation_type == "transcript":
            simplify = False
    if not annotation:
        simplify = False

    num_row = 2 if dendro else 1
    num_col = 2 if annotation else 1

    subplots = plots.create_subplots(num_col, num_row)
    return window, dendro, annotation, simplify, num_row, num_col, subplots


def browser_information(
    args,
    button_confirm,
    button_o3,
    button_o10,
    button_i3,
    button_i10,
    input_box,
    hierarchical_clustering,
    annotation,
    annotation_type,
    windowregion,
):
    if windowregion is None and input_box["props"]["value"] is None:
        return None, None, None, None, None, None, None
    window = window_input(input_box, windowregion)
    if window == (
        "error",
        "Input not recognized. Please enter genomic region in a valid format. Example chr20:58,839,718-58,911,192 or chr20:58839718-58911192",
    ):
        return None, None, None, None, None, None, None
    if window == ("The start position must be less than the end position."):
        return None, None, None, None, None, None, None
    if window == (
        "The region is too large. Please enter a region smaller than 1,000,000 bp."
    ):
        return None, None, None, None, None, None, None
    dendro = "on" in hierarchical_clustering
    annotation = "on" in annotation
    if annotation:
        if annotation_type == "gene":
            simplify = True
        if annotation_type == "transcript":
            simplify = False
    if not annotation:
        simplify = False

    num_row = 2 if dendro else 1
    num_col = 2 if annotation else 1

    subplots = plots.create_subplots(num_col, num_row)
    return window, dendro, annotation, simplify, num_row, num_col, subplots


def mod_freq_data(args, window, upload_data=None, filename=None, last_modified=None):
    return read_mods(args, window, upload_data, filename, last_modified)


def meth_browser(app, args, gff_file, annotation_dir):
    @app.callback(
        [
            Output(component_id="plot", component_property="children"),
            Output(component_id="error-message", component_property="children"),
        ],
        [
            Input(component_id="confirm-button", component_property="n_clicks"),
            Input(component_id="button-o3", component_property="n_clicks"),
            Input(component_id="button-o10", component_property="n_clicks"),
            Input(component_id="button-i3", component_property="n_clicks"),
            Input(component_id="button-i10", component_property="n_clicks"),
            Input(component_id="hierarchical_clustering", component_property="value"),
            Input(component_id="annotation", component_property="value"),
            Input(component_id="annotation-file", component_property="value"),
            Input(component_id="annotation-type", component_property="value"),
            Input(component_id="color_scale", component_property="value"),
            Input(component_id="intermediate-data", component_property="data"),
        ],
        [
            State(component_id="input-box", component_property="children"),
        ],
    )
    def update_meth_browser(
        button_confirm,
        button_o3,
        button_o10,
        button_i3,
        button_i10,
        hierarchical_clustering,
        annotation,
        annotation_file,
        annotation_type,
        color_scale,
        mod_data,
        input_box,
    ):
        if annotation == "off":
            annotation = "off"
            gff = None
        if annotation == "on":
            if annotation_file is None:
                annotation = "off"
                gff = None
            else:
                annotation = "on"
                gff = annotation_dir + "/" + annotation_file + "_sorted.gff3.gz"

        if input_box["props"]["value"] is None and mod_data is None:
            return None, None
        else:
            window, dendro, annotation, simplify, num_row, num_col, subplots = (
                browser_information(
                    args,
                    button_confirm,
                    button_o3,
                    button_o10,
                    button_i3,
                    button_i10,
                    input_box,
                    hierarchical_clustering,
                    annotation,
                    annotation_type,
                    args.window,
                )
            )
        if num_row is None and num_col is None:
            return None, None
        else:
            json_data = json.loads(mod_data)
            mod_data = pd.DataFrame(json_data["data"], columns=json_data["columns"])
            mod_data.index = json_data["index"]
            if mod_data.empty:
                return None, "No data for this region"
            # if args.gff:
            #     gff = args.gff
            # else:
            #     gff = gff_file
            fig = process_fig(
                mod_data,
                dendro,
                subplots,
                num_row,
                num_col,
                window,
                annotation,
                simplify,
                gff,
                color_scale,
                args.output,
            )
            return html.Div(dcc.Graph(figure=fig), id="plot"), None

    return html.Div(id="plot")


def dcc_store(app, args):

    @app.callback(
        [
            Output(component_id="intermediate-data", component_property="data"),
            Output(
                component_id="error-message-uploadowndata",
                component_property="children",
            ),
        ],
        [
            Input(component_id="confirm-button", component_property="n_clicks"),
            Input(component_id="button-o3", component_property="n_clicks"),
            Input(component_id="button-o10", component_property="n_clicks"),
            Input(component_id="button-i3", component_property="n_clicks"),
            Input(component_id="button-i10", component_property="n_clicks"),
            Input(component_id="upload-data", component_property="contents"),
            Input(component_id="input-box", component_property="children"),
        ],
        [
            State(component_id="input-box", component_property="children"),
            State(component_id="upload-data", component_property="filename"),
            State(component_id="upload-data", component_property="last_modified"),
        ],
    )
    def generate_data(
        confirm_button,
        button_o3,
        button_o10,
        button_i3,
        button_i10,
        upload_data,
        input_box,
        input_box_2,
        filename,
        last_modified,
    ):
        if (
            input_box["props"]["value"] is None
            and args.table is None
            and args.files is None
            and upload_data is None
        ):
            return None, None
        if input_box["props"]["value"] is None and upload_data is not None:
            args_False = False
            window = None
            mod_data = mod_freq_data(
                args_False, window, upload_data, filename, last_modified
            )
            mod_data.drop(columns=["chrom"], inplace=True)
            mod_data.set_index("position", inplace=True)
            json_data = mod_data.to_json(orient="split")
            return json_data, None
        else:
            window = window_input(input_box, args.window)
            if window == (
                "error",
                "Input not recognized. Please enter genomic region in a valid format. Example chr20:58,839,718-58,911,192 or chr20:58839718-58911192",
            ):
                return (
                    None,
                    "Input not recognized. Please enter genomic region in a valid format. Example chr20:58,839,718-58,911,192 or chr20:58839718-58911192",
                )
            if window == ("The start position must be less than the end position."):
                return None, "The start position must be less than the end position."
            if window == (
                "The region is too large. Please enter a region smaller than 1,000,000 bp."
            ):
                return (
                    None,
                    "The region is too large. Please enter a region smaller than 1,000,000 bp.",
                )
            window_region = Region(window)
            if upload_data is None:
                mod_data = mod_freq_data(args, window_region)
                mod_data = mod_data.reset_index(level="chrom", drop=True)
                json_data = mod_data.to_json(orient="split")
            elif upload_data is not None:
                args_False = False
                mod_data = mod_freq_data(
                    args_False, window_region, upload_data, filename, last_modified
                )
                if mod_data.empty:
                    return None, "No data for this region"
                else:
                    json_data = mod_data.to_json(orient="split")
            return json_data, None

    return html.Div([dcc.Store(id="intermediate-data")])


def window_input(
    input_box,
    window,
):
    if input_box["props"]["value"] is None and window is None:
        window = None
    else:
        window_input = input_box["props"]["value"] if input_box else window
        if (
            validate_input(window_input)
            == "The start position must be less than the end position."
        ):
            return "The start position must be less than the end position."
        if (
            validate_input(window_input)
            == "The region is too large. Please enter a region smaller than 1,000,000 bp."
        ):
            return "The region is too large. Please enter a region smaller than 1,000,000 bp."
        if not validate_input(window_input):
            return (
                "error",
                "Input not recognized. Please enter genomic region in a valid format. Example chr20:58,839,718-58,911,192 or chr20:58839718-58911192",
            )
    return window_input


def window_input_1000Genomes(args, input_box, genes_to_coords, window):
    if input_box["props"]["value"] is None and window is None:
        window = None
    else:
        window_input = input_box["props"]["value"] if input_box else window
        if (
            validate_input_1000Genomes(window_input, args)
            == "The start position must be less than the end position."
        ):
            return "The start position must be less than the end position."
        if (
            validate_input_1000Genomes(window_input, args)
            == "The region is too large. Please enter a region smaller than 1,000,000 bp."
        ):
            return "The region is too large. Please enter a region smaller than 1,000,000 bp."
        if (
            validate_input_1000Genomes(window_input, args)
            == "Invalid chromosome. Chromsome not recognized/not present in the data."
        ):
            return (
                "Invalid chromosome. Chromsome not recognized/not present in the data."
            )
        if (
            validate_input_1000Genomes(window_input, args)
            == "No data found for the given interval."
        ):
            return "No data found for the given interval."
        if not validate_input_1000Genomes(window_input, args):
            window_input = window_input.upper()
            coords = genes_to_coords.get(window_input)

            if not coords:
                return (
                    "error",
                    "Input not recognized. Please enter genomic region in a valid format. Example chr20:58,839,718-58,911,192 or chr20:58839718-58911192",
                )
            else:
                window_input = coords

    return window_input


def validate_input(input_text):
    input_text = input_text.replace(",", "")
    if ":" not in input_text or "-" not in input_text:
        return False
    chrom, positions = input_text.split(":")
    start, end = map(int, positions.split("-"))
    difference = end - start
    if difference < 0:
        return "The start position must be less than the end position."
    if difference > 1000000:
        return (
            "The region is too large. Please enter a region smaller than 1,000,000 bp."
        )
    return True


def validate_input_1000Genomes(input_text, args):
    db = args.db
    tbi_file = db + ".tbi"
    if not os.path.exists(tbi_file):
        sys.exit(f"ERROR: {tbi_file} 1000Genomes --db input not found")
    tabix_file = pysam.TabixFile(db)
    input_text = input_text.replace(",", "")
    if ":" not in input_text or "-" not in input_text:
        return False
    chrom, positions = input_text.split(":")
    start, end = map(int, positions.split("-"))
    difference = end - start
    if difference < 0:
        return "The start position must be less than the end position."
    if difference > 1000000:
        return (
            "The region is too large. Please enter a region smaller than 1,000,000 bp."
        )

    valid_chromosomes = set()
    for chromosome in tabix_file.contigs:
        valid_chromosomes.add(chromosome)
    if chrom not in valid_chromosomes:
        return "Invalid chromosome. Chromsome not recognized/not present in the data."
    try:
        records = tabix_file.fetch(chrom, start, end)
        for record in records:
            return True
        return "No data found for the given interval."
    except ValueError:
        return False


def make_gene_to_coords_dict(gff_file):
    genes_to_coords = {}
    if gff_file.endswith(".gz"):
        open_func = gzip.open
        mode = "rt"
    else:
        open_func = open
        mode = "r"
    with open_func(gff_file, mode) as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields[2] == "gene":
                name = [f for f in fields[8].split(";") if f.startswith("gene_name")][
                    0
                ].split("=")[1]
                genes_to_coords[name] = f"{fields[0]}:{fields[3]}-{fields[4]}"
    return genes_to_coords


def input_box(app, args):
    @app.callback(
        Output(
            component_id="input-box",
            component_property="children",
        ),
        Output(
            component_id="error-message-inputbox",
            component_property="children",
        ),
        [
            Input(component_id="confirm-button", component_property="n_clicks"),
            Input(component_id="button-o3", component_property="n_clicks"),
            Input(component_id="button-o10", component_property="n_clicks"),
            Input(component_id="button-i3", component_property="n_clicks"),
            Input(component_id="button-i10", component_property="n_clicks"),
        ],
        [
            State(component_id="input-box", component_property="children"),
        ],
    )
    def update_value(
        confirm_button, button_o3, button_o10, button_i3, button_i10, window
    ):
        if window["props"]["value"] is not None or args.window is not None:
            window = window["props"]["value"] if window else args.window
            if not validate_input(window):
                window = "error"
                return (
                    html.Div(
                        dcc.Input(type="text", value=window),
                        id="input-box",
                        style={"height": "30px", "margin": "0px 2px"},
                    ),
                    None,
                )
            window = Region(window)
            if "button-o3" == ctx.triggered_id:
                window = window * 3
            elif "button-o10" == ctx.triggered_id:
                window = window * 10
            elif "button-i3" == ctx.triggered_id:
                window = window / 3
            elif "button-i10" == ctx.triggered_id:
                window = window / 10
            window = window.fmt
            return (
                html.Div(
                    dcc.Input(type="text", value=window),
                    id="input-box",
                    style={
                        "height": "30px",
                        "margin": "0px 2px",
                    },
                ),
                None,
            )

        elif window["props"]["value"] is None and args.window is None:
            window = None
            return (
                html.Div(
                    dcc.Input(type="text", value=window),
                    id="input-box",
                    style={"height": "30px"},
                ),
                None,
            )

    return html.Div(
        dcc.Input(type="text", value=args.window),
        id="input-box",
        style={
            "height": "30px",
            "margin": "0px 2px",
        },
    )


if __name__ == "__main__":
    main()
