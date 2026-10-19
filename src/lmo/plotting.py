"""Shared matplotlib styling: one colour per C-rate, 300 dpi PNG output."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from . import FIGURES_DIR, RATES  # noqa: E402

COLORS = {rate: meta["color"] for rate, meta in RATES.items()}
CV_COLORS = {"0.1 mV/s": "#1f77b4", "0.2 mV/s": "#d62728", "0.3 mV/s": "#2ca02c"}
MEASURED_STYLE = dict(linestyle="-", linewidth=1.8, alpha=0.9)
MODEL_STYLE = dict(linestyle="--", linewidth=1.6)

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "figure.dpi": 100,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "legend.frameon": False,
    }
)


def savefig(fig, name: str):
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    path = FIGURES_DIR / f"{name}.png"
    fig.savefig(path)
    plt.close(fig)
    return path
