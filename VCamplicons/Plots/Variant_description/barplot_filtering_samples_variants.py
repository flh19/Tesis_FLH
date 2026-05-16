#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===========================================================================
Filtering Summary Figure
===========================================================================

This script generates publication-quality stacked barplots summarizing
sample and variant filtering across populations/cohorts.

The figure contains two panels:

(a) Percentage of samples retained and removed after filtering
(b) Percentage of variants retained and removed after filtering

INPUT FILE
----------
Tab-separated file containing:

Poblacion              Cohort name
muestras_antes         Number of samples before filtering
muestras_despues       Number of samples after filtering
variantes_antes        Number of variants before filtering
variantes_despues      Number of variants after filtering

Example:

Poblacion  muestras_antes muestras_despues variantes_antes variantes_despues
Di@bet.es  4618           4435             296              247
Hortega    1301           1265             119              95

FEATURES
--------
- Publication-quality figures
- Automatic normalization to percentages
- PNG and PDF outputs
- Configurable aesthetics
- External configuration file (no hardcoded counts)

OUTPUTS
-------
figuras_filtros/
├── figura_filtros_muestras_variantes.png
└── figura_filtros_muestras_variantes.pdf

USAGE
-----
python3 filtering_summary.py \
    --infile filtros_muestras_variantes.txt

=========================================================================== 
"""

# ==========================================================================
# IMPORTS
# ==========================================================================

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.ticker import FuncFormatter


# ==========================================================================
# GLOBAL CONFIGURATION
# ==========================================================================

OUTPUT_DIR = Path("figuras_filtros")
OUTPUT_DIR.mkdir(exist_ok=True)

DPI = 600

BAR_WIDTH = 0.65

EDGE_COLOR = "black"

COLOR_BEFORE = "#B0BEC5"
COLOR_AFTER = "#4682B4"

plt.rcParams.update({

    "font.family": "sans-serif",
    "axes.linewidth": 1.0
})


# ==========================================================================
# PERCENTAGE FORMATTER
# ==========================================================================

def percentage_formatter(x, pos):
    """
    Format y-axis labels as percentages.
    """

    return f"{round(x)}"


y_formatter = FuncFormatter(percentage_formatter)


# ==========================================================================
# AUXILIARY FUNCTIONS
# ==========================================================================

def normalize_to_percentage(values_after, values_before):
    """
    Convert absolute counts into percentages.
    """

    return (
        np.array(values_after) /
        np.array(values_before)
    ) * 100


def style_axes(ax, ymax=None):
    """
    Apply standardized styling to axes.
    """

    if ymax is None:
        ax.set_ylim(0, 100)

    else:
        ax.set_ylim(0, ymax * 1.05)

    ax.yaxis.set_major_formatter(y_formatter)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def adjust_fonts(
    ax,
    label_size=18,
    tick_size=18,
    legend_size=18
):
    """
    Standardize font sizes.
    """

    ax.xaxis.label.set_size(label_size)
    ax.yaxis.label.set_size(label_size)

    ax.tick_params(
        axis="both",
        labelsize=tick_size
    )

    legend = ax.get_legend()

    if legend is not None:

        for text in legend.get_texts():
            text.set_fontsize(legend_size)


def save_figure(fig, filename):
    """
    Save figure in PNG and PDF formats.
    """

    fig.savefig(
        OUTPUT_DIR / f"{filename}.png",
        dpi=DPI
    )

    fig.savefig(
        OUTPUT_DIR / f"{filename}.pdf",
        dpi=DPI
    )

    plt.close(fig)


# ==========================================================================
# MAIN FUNCTION
# ==========================================================================

def main():

    # ----------------------------------------------------------------------
    # ARGUMENT PARSER
    # ----------------------------------------------------------------------

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--infile",
        required=True,
        help="Input TSV file containing filtering statistics."
    )

    args = parser.parse_args()

    # ----------------------------------------------------------------------
    # LOAD INPUT DATA
    # ----------------------------------------------------------------------

    df = pd.read_csv(
        args.infile,
        sep="\t"
    )

    # ----------------------------------------------------------------------
    # EXTRACT VARIABLES
    # ----------------------------------------------------------------------

    poblaciones = df["Poblacion"].to_numpy()

    muestras_antes = df["muestras_antes"].to_numpy()
    muestras_despues = df["muestras_despues"].to_numpy()

    variantes_antes = df["variantes_antes"].to_numpy()
    variantes_despues = df["variantes_despues"].to_numpy()

    # ----------------------------------------------------------------------
    # NORMALIZATION
    # ----------------------------------------------------------------------
    # Percentages are calculated relative to pre-filtering values

    muestras_despues_pct = normalize_to_percentage(
        muestras_despues,
        muestras_antes
    )

    muestras_eliminadas_pct = 100 - muestras_despues_pct

    variantes_despues_pct = normalize_to_percentage(
        variantes_despues,
        variantes_antes
    )

    variantes_eliminadas_pct = 100 - variantes_despues_pct

    # ----------------------------------------------------------------------
    # X POSITIONS
    # ----------------------------------------------------------------------

    x = np.arange(len(poblaciones))

    # ----------------------------------------------------------------------
    # FIGURE INITIALIZATION
    # ----------------------------------------------------------------------

    fig, (ax1, ax2) = plt.subplots(
        1,
        2,
        figsize=(10, 4)
    )

    # ======================================================================
    # PANEL (A): MUESTRAS
    # ======================================================================

    ax1.bar(
        x,
        muestras_despues_pct,
        BAR_WIDTH,

        color=COLOR_AFTER,
        edgecolor=EDGE_COLOR,

        label="Posfiltrado"
    )

    ax1.bar(
        x,
        muestras_eliminadas_pct,
        BAR_WIDTH,

        bottom=muestras_despues_pct,

        color=COLOR_BEFORE,
        edgecolor=EDGE_COLOR,

        label="Prefiltrado"
    )

    ax1.set_ylabel("Muestras (%)")

    ax1.set_xticks(x)
    ax1.set_xticklabels(poblaciones)

    style_axes(ax1)

    adjust_fonts(ax1)

    ax1.text(
        -0.23,
        1.16,
        "(a)",

        transform=ax1.transAxes,

        fontsize=14,
        fontweight="bold",

        va="top",
        ha="left"
    )

    # ======================================================================
    # PANEL (B): VARIANTES
    # ======================================================================

    ax2.bar(
        x,
        variantes_despues_pct,
        BAR_WIDTH,

        color=COLOR_AFTER,
        edgecolor=EDGE_COLOR,

        label="Posfiltrado"
    )

    ax2.bar(
        x,
        variantes_eliminadas_pct,
        BAR_WIDTH,

        bottom=variantes_despues_pct,

        color=COLOR_BEFORE,
        edgecolor=EDGE_COLOR,

        label="Prefiltrado"
    )

    ax2.set_ylabel("Variantes (%)")

    ax2.set_xticks(x)
    ax2.set_xticklabels(poblaciones)

    style_axes(ax2)

    adjust_fonts(ax2)

    ax2.text(
        -0.23,
        1.16,
        "(b)",

        transform=ax2.transAxes,

        fontsize=14,
        fontweight="bold",

        va="top",
        ha="left"
    )

    # ======================================================================
    # SHARED LEGEND
    # ======================================================================

    handles, labels = ax1.get_legend_handles_labels()

    handles = handles[::-1]
    labels = labels[::-1]

    fig.legend(
        handles,
        labels,

        loc="lower center",

        bbox_to_anchor=(0.5, 0.001),

        ncol=2,

        frameon=False,

        fontsize=16
    )

    # ----------------------------------------------------------------------
    # FINAL LAYOUT
    # ----------------------------------------------------------------------

    fig.tight_layout(
        rect=[0, 0.10, 1, 1]
    )

    fig.subplots_adjust(
        wspace=0.26
    )

    # ----------------------------------------------------------------------
    # SAVE FIGURE
    # ----------------------------------------------------------------------

    save_figure(
        fig,
        "figura_filtros_muestras_variantes"
    )

    # ----------------------------------------------------------------------
    # FINAL REPORT
    # ----------------------------------------------------------------------

    print("\n====================================================")
    print("Filtering summary figure successfully generated")
    print("====================================================")

    print("\nOutput files:")

    print(
        f"  {OUTPUT_DIR / 'figura_filtros_muestras_variantes.png'}"
    )

    print(
        f"  {OUTPUT_DIR / 'figura_filtros_muestras_variantes.pdf'}"
    )


# ==========================================================================
# SCRIPT EXECUTION
# ==========================================================================

if __name__ == "__main__":
    main()
