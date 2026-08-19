"""Frozen figure style shared by every notebook in the repository.

This module centralizes the typography, palette, and line/marker defaults of
the manuscript figures. Notebooks call :func:`apply_paper_style` (main text) or
:func:`apply_si_style` (Supporting Information) and never set these values
locally, so that the two sets of figures stay visually consistent.
"""

from __future__ import annotations

from cycler import cycler
import matplotlib as mpl

# Increment this label only through an intentional paper-wide style decision.
MAIN_FIGURE_STYLE_VERSION = "pnas-main-v1"

PAPER_COLORS = ["#1A4E8A", "#B22222", "#E69F00"]
DIST_COLORS = {
    "patents": "#1A4E8A",
    "papers": "#E69F00",
    "sim": "#B22222",
    "reference": "#222222",
}

# Semantic palette used by the main empirical/comparison figures.
# Keep the mapping fixed when building SI diagnostics so that visual meaning
# does not change between the main text and supporting information.
OBSERVABLE_COLORS = {
    "novelty": "#1A4E8A",
    "intrinsic_time": "#B22222",
    "explorers": "#E69F00",
    "reference": "#222222",
}

POINT_SIZE = 10
LINE_WIDTH = 5

# Stroke weights of the SI panels. The SI figures use a wider canvas than the
# main ones, so they set line widths per artist instead of relying on
# ``lines.linewidth``; these are the three weights they use.
SI_DATA_LW = 3.0     # empirical / simulated curves
SI_FIT_LW = 3.0      # fitted curves drawn on top of the data
SI_REF_LW = 2.4      # reference lines and guides
SI_LEGEND_SIZE = 25  # SI legends run slightly larger than the main-text ones


def apply_paper_style(grid: bool = False) -> dict[str, object]:
    """Apply the paper-style matplotlib configuration.

    Parameters
    ----------
    grid:
        If True, enable a light background grid. This is useful for the
        diagnostic panels while preserving the paper palette and typography.
    """

    mpl.rcParams.update(
        {
            "font.size": 28,
            "axes.labelsize": 30,
            "axes.titlesize": 30,
            "xtick.labelsize": 28,
            "ytick.labelsize": 28,
            "legend.fontsize": 22,
            "figure.figsize": (10, 6),
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "axes.prop_cycle": cycler(color=PAPER_COLORS),
            "lines.linewidth": LINE_WIDTH,
            "lines.markersize": POINT_SIZE,
            "axes.linewidth": 1.8,
            "xtick.major.width": 1.6,
            "ytick.major.width": 1.6,
            "xtick.major.size": 6.0,
            "ytick.major.size": 6.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": grid,
            "grid.alpha": 0.2,
            "grid.linewidth": 1.0,
            "grid.color": "#C9CDD3",
            "grid.linestyle": "-",
            "agg.path.chunksize": 10000,
        }
    )
    return {
        "point_size": POINT_SIZE,
        "line_width": LINE_WIDTH,
        "series_colors": tuple(PAPER_COLORS),
        "dist_colors": dict(DIST_COLORS),
        "observable_colors": dict(OBSERVABLE_COLORS),
        "style_version": MAIN_FIGURE_STYLE_VERSION,
    }


def apply_si_style() -> dict[str, object]:
    """Apply the paper style with the SI legend size.

    Same typography and palette as the main figures; only the legend runs
    larger, because the SI panels are wider.
    """

    style = apply_paper_style(grid=False)
    mpl.rcParams["legend.fontsize"] = SI_LEGEND_SIZE
    return style
