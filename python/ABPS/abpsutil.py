import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os


def _3d_plot(
    x, y, z,
    xtitle_text, ytitle_text, ztitle_text,
    xrange=None, yrange=None, zrange=None,
    showscale=False, cmin=None, cmax=None,
    title="",
):
    if cmin is None:
        cmin = z.min()
    if cmax is None:
        cmax = z.max()
    if xrange is None:
        xrange = [x.min(), x.max()]
    if yrange is None:
        yrange = [y.min(), y.max()]
    if zrange is None:
        zrange = [cmin, cmax]
    fig = make_subplots(rows=1, cols=1, specs=[[{'is_3d': True}]],)
    fig.add_trace(go.Surface(z=z, x=x, y=y,
                             showscale=showscale, cmin=cmin, cmax=cmax), row=1, col=1)
    scene = dict(
        yaxis=dict(title_text=ytitle_text, range=yrange),
        xaxis=dict(title_text=xtitle_text, range=xrange),
        zaxis=dict(title_text=ztitle_text, range=zrange),
        camera=dict(
            eye=dict(
                x=1.5,
                y=1.5,
                z=1.5,
            )
        )
    )
    layout = go.Layout(title=title,title_x=0.5, scene=scene,
                       margin=dict(l=0, r=0, b=0, t=0, pad=0),
                       )
    fig.update(layout=layout)

    return fig


def _heatmap_plot(x, y, z, xtitle_text, ytitle_text, showscale=False,title=""):
    fig = make_subplots(rows=1, cols=1, specs=[[{'type': "heatmap"}]],)
    fig.add_trace(go.Heatmap(z=z, x=x, y=y, showscale=showscale,), row=1, col=1)
    layout = go.Layout(title=title,title_x=0.5,
        yaxis=dict(title=ytitle_text, title_text=ytitle_text),
        xaxis=dict(title=xtitle_text, title_text=xtitle_text),)
    fig.update(layout=layout)

    return fig


def save_img(fig, path, name):

    sm = dict(width=500, height=400, scale=3)
    #md = dict(width=600, height=600, scale=1)
    #lg = dict(width=800, height=500, scale=1)

    _fn = os.path.join(path, "sm", name)
    fig.write_image(_fn, **sm)

    #_fn = os.path.join(path, "md", name)
    #fig.write_image(_fn, **md)

    #_fn = os.path.join(path, "lg", name)
    #fig.write_image(_fn, **lg)

    # plotly.offline.plot(fig, filename = 'akhm_analytic_numeric.html', auto_open=False)
