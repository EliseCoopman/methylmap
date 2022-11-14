import plotly
import plotly.graph_objects as go


def make_per_gene_annot_line_trace(transcript, window, x_pos):
    """Generate a line trace for the gene
    Trace can get limited by the window sizes
    """
    return go.Scatter(
        y=[max(transcript.begin, window.begin), min(transcript.end, window.end)],
        x=[x_pos, x_pos],
        mode="lines",
        line=dict(width=2, color=transcript.color),
        name=transcript.transcript,
        text=transcript.gene,
        hoverinfo="text",
        showlegend=False,
    )


def make_per_exon_arrow_trace(transcript, begin, end, x_pos):
    """Generate a line+marker trace for the exon
    The shape is an arrow, as defined by the strand in transcript.marker
    """
    return go.Scatter(
        y=[begin,end],
        x=[x_pos, x_pos],
        mode="lines+markers",
        line=dict(width=8, color=transcript.color),
        name=transcript.transcript,
        text=transcript.gene,
        hoverinfo="text",
        showlegend=False,
        marker=dict(symbol=transcript.marker, size=8),
    )


def plot_methylation(subplots, meth_data, num_col):
    samplelist = list(meth_data)
    positionlist = meth_data.index.values.tolist()
    overviewarray = meth_data.to_numpy()

    fig = subplots.add_trace(
        go.Heatmap(z=overviewarray, x=samplelist, y=positionlist), row=1, col=num_col
    )
    fig.update_xaxes(tickangle=45, tickfont=dict(size=10), row=1, col=num_col)
    fig.update_yaxes(tickfont=dict(size=10), row=1, col=num_col, autorange="reversed")
    return fig


def create_subplots(num_col):
    """
    Prepare the panels for the subplots in case of annotation track.
    """
    if num_col > 1:
        fig = plotly.subplots.make_subplots(
            rows=1,
            cols=num_col,
            column_widths=[0.1, 0.9],
            horizontal_spacing=0.001,
            shared_yaxes=True,
        )
    else:
        fig = plotly.subplots.make_subplots(rows=1, cols=1)
    fig.update_layout(
        {
            "plot_bgcolor": "rgba(0,0,0,0)",
            "paper_bgcolor": "rgba(0,0,0,0)",
        },
        title="methylation frequency",
    )
    return fig