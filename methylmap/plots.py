import sys
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

def create_output_methylmap(fig,outfig,window):
    if outfig is None:
        outfig = f'methylmap_{window.fmt}.html'
    else:
        from pathlib import Path
        outfig = outfig.format(region=window.fmt)
        p = Path(outfig)
        Path.mkdir(p.parent, exist_ok=True,parents=True)

    create_output(fig, outfig)

def create_output_dendrogram(fig, outdendro, window):
    if outdendro is None:
        outdendro = f'dendrogram_{window.fmt}.html'
    else:
        from pathlib import Path
        outdendro = outdendro.format(region=window.fmt)
        p = Path(outdendro)
        Path.mkdir(p.parent, exist_ok=True,parents=True)

    create_output(fig, outdendro)

def create_output(fig,outfig):
    if outfig.endswith(".html"):
        html_output(fig,outfig)
    else:
        try:
            fig.write_image(outfig)
        except ValueError as e:
            sys.stderr.write("\n\nERROR: creating the image in this file format failed.\n")
            sys.stderr.write("ERROR: creating in default html format instead.\n")
            sys.stderr.write("ERROR: additional packages required. Detailed error:\n")
            sys.stderr.write(str(e))
            html_output(fig, outfig)
    
def html_output(fig, outfig):
    with open(outfig, "w+") as output:
        output.write(plotly.offline.plot(fig, output_type="div", show_link=False,include_plotlyjs="cdn"))