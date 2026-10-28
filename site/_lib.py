"""Shared helpers for the website: data loading, identity colors, Plotly theme."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio

ROOT = Path.cwd()
if not (ROOT / "config").exists():  # rendered with execute-dir: project (site/)
    ROOT = ROOT.parent
FITSET = os.environ.get("LEAFRTM_FITSET", "pilot")
# Site data: results/<fitset>/site_data (built by the workflow, or unpacked from
# the GitHub release asset site-data-<fitset>.tar.zst).
DATA = Path(os.environ.get("LEAFRTM_SITE_DATA", ROOT / "results" / FITSET / "site_data"))

# ---- Identity: one fixed color per RTM, validated categorical order (light mode);
# marker symbols are the secondary encoding.
MODELS = ["prospect4", "prospect5", "prospect5b", "prospectd", "prospectpro", "liberty", "leafsip"]
MODEL_LABEL = {"prospect4": "PROSPECT-4", "prospect5": "PROSPECT-5", "prospect5b": "PROSPECT-5B",
               "prospectd": "PROSPECT-D", "prospectpro": "PROSPECT-PRO", "liberty": "LIBERTY (Eq. 16 fixed)",
               "leafsip": "LeafSIP"}
MODEL_COLOR = dict(zip(MODELS, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
                                "#008300", "#6250d6"]))
SYMBOLS = ["circle", "square", "diamond", "triangle-up", "x", "cross", "star"]
MODEL_SYMBOL = dict(zip(MODELS, SYMBOLS))
# LIBERTY as published (Dawson et al. 1998, Eq. 16 bug included): a reference, not fitted;
# muted gray so it never competes with the fitted models' hues.
MODEL_LABEL["liberty_dawson"] = "LIBERTY (as published)"
MODEL_COLOR["liberty_dawson"] = "#898781"
MODEL_SYMBOL["liberty_dawson"] = "cross-open"
ERR_LABEL = {"ou": "Hetero + OU", "ss": "Hetero + Matérn-3/2 (fixed ℓ)",
             "ld": "Learned discrepancy"}
ERR_COLOR = {"ou": "#2a78d6", "ss": "#eb6834", "ld": "#1baf7a"}
STATE_COLOR = {"fresh": "#2a78d6", "dried_intact": "#eb6834", "pressed": "#1baf7a",
               "ground": "#eda100"}
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE = "#e8e7e2", "#c3c2b7", "#fcfcfb"

TRAIT = {  # label, unit
    "chl_area": ("Chlorophyll a+b", "µg cm⁻²"), "car_area": ("Carotenoids", "µg cm⁻²"),
    "anth_area": ("Anthocyanins", "µg cm⁻²"), "ewt": ("Water (EWT)", "g m⁻²"),
    "lma": ("Leaf mass per area", "g m⁻²"), "n_area": ("Nitrogen per area", "g m⁻²"),
    "protein_area": ("Protein per area", "g m⁻²"),
    "lignin_cellulose_area": ("Lignin + cellulose per area", "g m⁻²"),
}


def trait_label(t, unit=True):
    if t in TRAIT:
        return f"{TRAIT[t][0]} ({TRAIT[t][1]})" if unit else TRAIT[t][0]
    return t.replace("_", " ")


def run_label(r) -> str:
    """'PROSPECT-D · Learned discrepancy · R/T-correlated · R only'."""
    parts = [MODEL_LABEL.get(r["model"], r["model"]), ERR_LABEL.get(r["errmodel"], r["errmodel"])]
    if bool(r.get("rt")):
        parts.append("R/T-correlated")
    if r.get("data_use") not in (None, "all"):
        parts.append(f"{r['data_use']} only")
    return " · ".join(parts)


pio.templates["leafrtm"] = go.layout.Template(layout=dict(
    font=dict(family='system-ui, -apple-system, "Segoe UI", sans-serif', color=INK, size=13),
    paper_bgcolor=SURFACE, plot_bgcolor=SURFACE,
    xaxis=dict(gridcolor=GRID, linecolor=AXIS, zerolinecolor=AXIS, tickfont=dict(color=INK2),
               title=dict(font=dict(color=INK2))),
    yaxis=dict(gridcolor=GRID, linecolor=AXIS, zerolinecolor=AXIS, tickfont=dict(color=INK2),
               title=dict(font=dict(color=INK2))),
    legend=dict(font=dict(color=INK2), bgcolor="rgba(0,0,0,0)"),
    hoverlabel=dict(bgcolor="white", font=dict(color=INK)),
    margin=dict(l=60, r=20, t=50, b=50), hovermode="closest",
))
pio.templates.default = "leafrtm"


def show(fig, height=None):
    """Render a Plotly figure using the page's local plotly.min.js (offline)."""
    from IPython.display import HTML, display
    if height:
        fig.update_layout(height=height)
    lg = fig.layout.legend
    if lg.orientation == "h" and (lg.y or 0) > 1:  # top legend: pin to the container, above subplot titles
        fig.update_layout(legend=dict(yref="container", y=0.99, yanchor="top", x=0, xanchor="left"),
                          margin=dict(t=max(fig.layout.margin.t or 0, 110)))
    if fig.layout.sliders:  # keep slider labels clear of the x-axis titles
        fig.update_layout(margin=dict(b=150))
        for sl in fig.layout.sliders:
            sl.pad = dict(t=75)
    display(HTML(fig.to_html(full_html=False, include_plotlyjs=False,
                             config={"displaylogo": False, "responsive": True})))


try:  # itables bundles its JS into the page
    from itables import init_notebook_mode
    init_notebook_mode(all_interactive=False, connected=False)
except Exception:  # pragma: no cover
    pass


def show_table(df, **kw):
    try:
        from itables import show as ishow
        opts = {"classes": "display compact", "maxBytes": 0, "scrollX": True, "showIndex": False, **kw}
        ishow(df, **opts)
    except Exception:
        from IPython.display import display
        display(df)


def read(name: str) -> pd.DataFrame:
    p = DATA / f"{name}.parquet"
    return pd.read_parquet(p) if p.exists() else pd.DataFrame()


def fmt(df, digits=3):
    out = df.copy()
    for c in out.select_dtypes("number"):
        out[c] = out[c].map(lambda v: "" if pd.isna(v) else f"{v:.{digits}g}")
    return out


def stats(items):
    """Row of hero numbers: [(value, label), ...]."""
    from IPython.display import HTML, display
    cells = "".join(f'<div class="stat"><div class="v">{v}</div><div class="l">{lab}</div></div>'
                    for v, lab in items)
    display(HTML(f'<div class="stat-row">{cells}</div>'))


def log_ticks(lo, hi):
    lo, hi = max(lo, 1e-12), max(hi, 1e-12)
    vals = []
    for e in range(int(np.floor(np.log10(lo))), int(np.ceil(np.log10(hi))) + 1):
        for m in (1, 2, 5):
            v = m * 10.0**e
            if lo <= v <= hi:
                vals.append(v)
    return vals


def log_axis(fig, axis, values, row=None, col=None, title=None):
    """Log axis with readable 1-2-5 ticks spanning ``values``."""
    v = np.asarray(values, float)
    v = v[np.isfinite(v) & (v > 0)]
    if len(v) == 0:
        return
    ticks = log_ticks(v.min(), v.max())
    if len(ticks) < 3:  # narrow range: add intermediate ticks
        ticks = sorted(set(ticks) | {float(f"{x:.2g}") for x in np.geomspace(v.min(), v.max(), 5)})
    kw = dict(type="log", tickvals=ticks, ticktext=[f"{t:g}" for t in ticks])
    if title:
        kw["title_text"] = title
    upd = fig.update_xaxes if axis == "x" else fig.update_yaxes
    upd(**kw, **({"row": row, "col": col} if row else {}))


def one_to_one(fig, lo, hi, row=None, col=None):
    kw = dict(row=row, col=col) if row else {}
    fig.add_trace(go.Scatter(x=[lo, hi], y=[lo, hi], mode="lines", line=dict(color=AXIS, width=1),
                             hoverinfo="skip", showlegend=False), **kw)


def no_data(msg="No results yet for this fit set."):
    from IPython.display import Markdown, display
    display(Markdown(f"::: {{.callout-note}}\n{msg}\n:::"))
