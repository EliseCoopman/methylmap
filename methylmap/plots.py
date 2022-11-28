import logging
import plotly.subplots
import plotly.graph_objects as go


def create_subplots(num_col):
    """
    Prepare the panels for the subplots in case of annotation track.
    """
    if num_col > 1:
        logging.info("When annotation input, make 2 subplots")
        fig = plotly.subplots.make_subplots(
            rows=1,
            cols=num_col,
            column_widths=[0.1, 0.9],
            horizontal_spacing=0.001,
            shared_yaxes=True
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


def plot_methylation(subplots, meth_data, num_col):
    """Make heatmap of modification frequencies."""

    samplelist = list(meth_data)
    positionlist = meth_data.index.values.tolist()
    overviewarray = meth_data.to_numpy()
    
    fig = subplots.add_trace(   
    go.Heatmap(z=overviewarray, x=samplelist, y=positionlist), row=1, col=num_col
    )
    
    fig.update_xaxes(tickangle=45, tickfont=dict(size=4), row=1, col=num_col)
    fig.update_yaxes(tickfont=dict(size=10), row=1, col=num_col, autorange="reversed")
    
    return fig
