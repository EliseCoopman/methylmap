import sys
import logging
import plotly.subplots
import plotly.graph_objects as go


def create_subplots(num_col, num_row):
    """
    Prepare the panels for the subplots in case of annotation track.
    """
    if num_col == 2 and num_row == 2:
        logging.info("When annotation input and/or dendro input, make subplots")
        fig = plotly.subplots.make_subplots(
            rows=num_row,
            cols=num_col,
            column_width=[0.1, 0.9],
            row_width=[0.7, 0.3],
            horizontal_spacing=0.00000,
            vertical_spacing=0.00000,
            shared_xaxes=True,
            shared_yaxes=True,
        )
    elif num_col == 2 and num_row == 1:
        fig = plotly.subplots.make_subplots(
            rows=num_row,
            cols=num_col,
            column_width=[0.1, 0.9],
            horizontal_spacing=0.00000,
            vertical_spacing=0.00000,
            shared_yaxes=True,
        )
    elif num_col == 1 and num_row == 2:
        fig = plotly.subplots.make_subplots(
            rows=num_row,
            cols=num_col,
            row_width=[0.7, 0.3],
            horizontal_spacing=0.00000,
            vertical_spacing=0.00000,
            shared_xaxes=True,
        )
    else:
        fig = plotly.subplots.make_subplots(rows=1, cols=1)
    fig.update_layout(
        {
            "plot_bgcolor": "rgba(0,0,0,0)",
            "paper_bgcolor": "rgba(0,0,0,0)",
        }
    )
    return fig


def plot_methylation(subplots, meth_data, num_col, num_row, color_scale):
    """Make heatmap of modification frequencies."""

    samplelist = list(meth_data)
    positionlist = meth_data.index.values.tolist()
    overviewarray = meth_data.to_numpy()

    fig = subplots.add_trace(
        go.Heatmap(z=overviewarray, x=samplelist, y=positionlist),
        row=num_row,
        col=num_col,
    )

    fig.update_yaxes(
        tickfont=dict(size=10), row=num_row, col=num_col, autorange="reversed"
    )
    fig.update_xaxes(tickangle=45, tickfont=dict(size=4), row=num_row, col=num_col)
    cmaps = {
        "plasma": "Plasma",
        "greys": "Greys",
        "bluered": "Bluered"
    }
    fig.update_traces(colorscale=cmaps.get(color_scale, "Plasma"))

    return fig


def create_output_methylmap(fig, outfig, window):
    if outfig is None:
        outfig = f"methylmap_{window.fmt}.html"
    else:
        from pathlib import Path

        outfig = outfig.format(region=window.fmt)
        p = Path(outfig)
        Path.mkdir(p.parent, exist_ok=True, parents=True)

    create_output(fig, outfig)


def create_output(fig, outfig):
    if outfig.endswith(".html"):
        html_output(fig, outfig)
    else:
        try:
            fig.write_image(outfig)
        except ValueError as e:
            sys.stderr.write(
                "\n\nERROR: creating the image in this file format failed.\n"
            )
            sys.stderr.write("ERROR: creating in default html format instead.\n")
            sys.stderr.write("ERROR: additional packages required. Detailed error:\n")
            sys.stderr.write(str(e))
            html_output(fig, outfig)


def html_output(fig, outfig):
    with open(outfig, "w+") as output:
        output.write(
            plotly.offline.plot(
                fig, output_type="div", show_link=False, include_plotlyjs="cdn"
            )
        )
