import methylmap.plots as plots
import methylmap.annotation as annot
import methylmap.dendro as dendrogram
from methylmap.import_data import read_mods
from methylmap.region import Region
from methylmap.version import __version__

from dash import Dash, dcc, html, Input, Output, State, ctx


import os
import sys
import re
import numpy as np
from argparse import ArgumentParser


def main():
    args = get_args()
    app = Dash(__name__)
    app.layout = html.Div(
        children=[
            html.H1(
                children="Methylmap",
            ),
            input_box(app, args),
            html.Button(id="confirm-button", n_clicks=0, children="Confirm"),
            html.Button(id="button-o3", n_clicks=0, children="Zoom out 3x"),
            html.Button(id="button-o10", n_clicks=0, children="Zoom out 10x"),
            html.Button(id="button-i3", n_clicks=0, children="Zoom in 3x"),
            html.Button(id="button-i10", n_clicks=0, children="Zoom in 10x"),
            html.Div(
                [  # Wrap Annotation and annotation-type in a common div
                    dcc.Checklist(
                        id="annotation",
                        options=[{"label": "Annotation", "value": "annotation"}],
                        value=[],
                        style={"display": "inline-block"},
                    ),
                    dcc.RadioItems(
                        id="annotation-type",
                        options=[
                            {"label": "Gene", "value": "gene"},
                            {"label": "Transcript", "value": "transcript"},
                        ],
                        value="gene",
                        labelStyle={"display": "inline-block"},
                        style={"display": "inline-block"},
                    ),
                ]
            ),
            dcc.Checklist(
                id="hierarchical_clustering",
                options=[{"label": "Hierarchical clustering", "value": "dendro"}],
                value=[],
            ),
            html.Div(id="error-message", style={"color": "red"}),
            meth_browser(app, args),
        ]
    )

    @app.callback(Output("annotation-type", "style"), Input("annotation", "value"))
    def toggle_annotation_type_visibility(annotation):
        if "annotation" in annotation:
            return {"display": "block"}  # Show RadioItems when 'Annotation' is checked
        else:
            return {"display": "none"}  # Hide RadioItems when 'Annotation' is unchecked

    app.run_server(debug=True)
    # app.run()


def get_args():
    parser = ArgumentParser(
        description="Plotting tool for population scale nucleotide modifications."
    )
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument(
        "-f",
        "--files",
        nargs="+",
        help="list with nanopolish (processed with calculate_methylation_frequency.py) files or BAM/CRAM files",
    )
    action.add_argument("-t", "--table", help="methfreqtable or overviewtable input")
    parser.add_argument(
        "-w",
        "--window",
        help="region to visualise, format: chr:start-end (example: chr20:58839718-58911192)",
    )
    parser.add_argument("-n", "--names", nargs="*", default=[], help="list with sample names")
    parser.add_argument("--gff", "--gtf", help="add annotation track based on GTF/GFF file")
    parser.add_argument(
        "--expand",
        help="number of base pairs to expand the window with in both directions",
        type=int,
        default=0,
    )
    parser.add_argument("--outtable", help="file to write the frequencies table to in tsv format")
    parser.add_argument(
        "--outfig",
        help="file to write output heatmap to, default: methylmap_{chr}_{start}_{end}.html (missing paths will be created)",
    )
    parser.add_argument("--groups", nargs="*", help="list of experimental group for each sample")
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
        help="modified base of interest when BAM/CRAM files as input. Options are: 5mC, 5hmC, 5fC, 5caC, 5hmU, 5fU, 5caU, 6mA, 5oxoG, Xao, default = 5mC",
        default="5mC",
        choices=["5mC", "5hmC", "5fC", "5caC", "5hmU", "5fU", "5caU", "6mA", "5oxoG", "Xao"],
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
    return args


def meth_browser(app, args):
    @app.callback(
        [
            Output(component_id="plot", component_property="children"),
            Output("error-message", "children"),
        ],
        [
            Input(component_id="confirm-button", component_property="n_clicks"),
            Input(component_id="button-o3", component_property="n_clicks"),
            Input(component_id="button-o10", component_property="n_clicks"),
            Input(component_id="button-i3", component_property="n_clicks"),
            Input(component_id="button-i10", component_property="n_clicks"),
            Input(component_id="hierarchical_clustering", component_property="value"),
            Input(component_id="annotation", component_property="value"),
            Input(component_id="annotation-type", component_property="value"),
        ],
        [State(component_id="input-box", component_property="children")],
    )
    def update_plots(
        button_confirm,
        button_o3,
        button_o10,
        button_i3,
        button_i10,
        hierarchical_clustering,
        annotation,
        annotation_type,
        window,
    ):
        window_input = window["props"]["value"] if window else args.window

        # Validate the input format
        if not validate_input(window_input):
            return (
                None,
                "Invalid input format. Please enter genomic region in a valid format. Example chr20:58,839,718-58,911,192 or chr20:58839718-58911192",
            )

        window = Region(window_input)
        if "button-o3" == ctx.triggered_id:
            window = window * 3
        elif "button-o10" == ctx.triggered_id:
            window = window * 10
        elif "button-i3" == ctx.triggered_id:
            window = window / 3
        elif "button-i10" == ctx.triggered_id:
            window = window / 10

        # Convert checklist values to boolean flags
        dendro = "dendro" in hierarchical_clustering
        annotation = "annotation" in annotation
        if annotation:
            if annotation_type == "gene":
                simplify = True
            if annotation_type == "transcript":
                simplify = False

        num_row = 2 if dendro else 1
        num_col = 2 if annotation else 1  # number of subplots (columns) needed

        subplots = plots.create_subplots(num_col, num_row)
        # frequencies table with all meth frequencies of all samples
        meth_data, window = read_mods(
            args.files,
            args.table,
            args.names,
            window,
            args.groups,
            args.gff,
            args.fasta,
            args.mod,
            args.hapl,
            dendro,
        )
        if dendro:
            meth_data, den, list_sorted_samples = dendrogram.make_dendro(meth_data, window)
        meth_data.to_csv(args.outtable, sep="\t", na_rep=np.NaN, header=True)
        fig = plots.plot_methylation(subplots, meth_data, num_col, num_row)

        if annotation:
            annotation_traces = annot.gff_annotation(args.gff, window, simplify)
            for annot_trace in annotation_traces:
                fig.append_trace(trace=annot_trace, row=num_row, col=1)
            fig.update_xaxes(
                title_text="", showticklabels=False, zeroline=False, row=num_row, col=1
            )
            fig.update_yaxes(title_text="", showticklabels=True, zeroline=False, row=num_row, col=1)

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
            # plots.create_output_methylmap(fig, outfig, window)
        return html.Div(dcc.Graph(figure=fig), id="plot"), None  # No error message

    return html.Div(id="plot")


def input_box(app, args):
    @app.callback(
        Output(component_id="input-box", component_property="children"),
        [
            Input(component_id="button-o3", component_property="n_clicks"),
            Input(component_id="button-o10", component_property="n_clicks"),
            Input(component_id="button-i3", component_property="n_clicks"),
            Input(component_id="button-i10", component_property="n_clicks"),
        ],
        [State(component_id="input-box", component_property="children")],
    )
    def update_value(button_o3, button_o10, button_i3, button_i10, window):
        window = Region(window["props"]["value"]) if window else Region(args.window)
        if "button-o3" == ctx.triggered_id:
            window = window * 3
        elif "button-o10" == ctx.triggered_id:
            window = window * 10
        elif "button-i3" == ctx.triggered_id:
            window = window / 3
        elif "button-i10" == ctx.triggered_id:
            window = window / 10

        return html.Div(dcc.Input(type="text", value=window.fmt), id="input-box")

    return html.Div(dcc.Input(type="text", value=args.window), id="input-box")


def validate_input(input_text):
    # Define a regular expression pattern for the valid formats
    pattern_with_commas = r"^chr\d+:\d+,\d+,\d+-\d+,\d+,\d+$"
    pattern_without_commas = r"^chr\d+:\d+-\d+$"

    # Use the re.match() function to check if the input matches either pattern
    if re.match(pattern_with_commas, input_text) or re.match(pattern_without_commas, input_text):
        return True
    else:
        return False


if __name__ == "__main__":
    main()
