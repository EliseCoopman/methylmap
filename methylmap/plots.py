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
        },
        title="Nucleotide modification frequencies",
    )
    return fig


def plot_methylation(meth_data, num_col, num_row):
    """Make heatmap of modification frequencies."""
    subplots = create_subplots(num_col, num_row)
    samplelist = list(meth_data)
    positionlist = meth_data.index.values.tolist()
    overviewarray = meth_data.to_numpy()

    fig = subplots.add_trace(
        go.Heatmap(z=overviewarray, x=samplelist, y=positionlist), row=num_row, col=num_col
    )

    fig.update_yaxes(tickfont=dict(size=10), row=num_row, col=num_col, autorange="reversed")
    fig.update_xaxes(tickangle=45, tickfont=dict(size=4), row=num_row, col=num_col)

    return fig
