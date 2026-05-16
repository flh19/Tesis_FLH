#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Per-amplicon mean coverage barplot for a target gene.

Description:
This script generates a publication-quality barplot showing
mean amplicon coverage across cohorts for a user-defined gene.

The workflow is fully generalizable to any gene by modifying
the GENE variable.
"""

# =========================================================
# Imports
# =========================================================

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt


# =========================================================
# Global plotting style
# =========================================================

plt.rcParams.update({

    "font.size": 10,

    "axes.labelsize": 15,
    "axes.titlesize": 12,

    "xtick.labelsize": 13,
    "ytick.labelsize": 15,

    "pdf.fonttype": 42,
    "ps.fonttype": 42
})


# =========================================================
# Global configuration
# =========================================================

# ---------------------------------------------------------
# Target gene
# ---------------------------------------------------------
# Modify this variable to analyze another gene

GENE = "ACAD10"

# ---------------------------------------------------------
# Input/output files
# ---------------------------------------------------------

INFILE = (
    f"nivel3_cobertura_media_por_amplicon_"
    f"{GENE}.tsv"
)

OUTPREFIX = (
    f"nivel3_barplot_cobertura_amplicon_"
    f"{GENE}"
)

# ---------------------------------------------------------
# Cohort color palette
# ---------------------------------------------------------

colores = {

    "Di@bet.es": "#4682B4",

    "Hortega": "#B0BEC5"
}


# =========================================================
# Data loading
# =========================================================

df = pd.read_csv(
    INFILE,
    sep="\t"
)

amplicones = sorted(
    df["amplicon"].unique()
)

poblaciones = list(
    colores.keys()
)

x = np.arange(
    len(amplicones)
)

width = 0.35


# =========================================================
# Figure initialization
# =========================================================

fig, ax = plt.subplots(
    figsize=(10, 4)
)


# =========================================================
# Barplot generation
# =========================================================

for i, pob in enumerate(poblaciones):

    # -----------------------------------------------------
    # Extract cohort-specific data
    # -----------------------------------------------------

    datos = (

        df[
            df["poblacion"] == pob
        ]

        .set_index("amplicon")

        .loc[amplicones]
    )

    # -----------------------------------------------------
    # Barplot
    # -----------------------------------------------------

    ax.bar(

        x - width / 2 + i * width,

        datos["cobertura_media"],

        width=width,

        color=colores[pob],

        alpha=0.8,

        edgecolor="black",

        linewidth=0.5,

        label=pob
    )


# =========================================================
# Axis formatting
# =========================================================

ax.set_xticks(x)

ax.set_xticklabels(
    "",
    rotation=0,
    ha="center"
)

ax.set_ylabel(
    "Mean coverage per base (X)"
)

ax.set_xlabel(
    GENE,
    fontstyle="italic"
)

ax.set_title("")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)


# =========================================================
# Optional legend
# =========================================================

# Uncomment to display legend

# handles = [
#     plt.Line2D(
#         [0],
#         [0],
#         color=colores[p],
#         lw=8,
#         label=p
#     )
#     for p in poblaciones
# ]
#
# ax.legend(
#     handles=handles,
#
#     frameon=False,
#
#     title="",
#
#     loc="upper right",
#
#     fontsize=14,
#
#     handlelength=2.5,
#
#     handleheight=1.5,
#
#     labelspacing=0.4,
#
#     borderpad=0.6
# )


# =========================================================
# Figure export
# =========================================================

plt.tight_layout()

plt.savefig(
    f"{OUTPREFIX}.png",
    dpi=600
)

plt.savefig(
    f"{OUTPREFIX}.pdf"
)

plt.close()


# =========================================================
# Final message
# =========================================================

print(
    "Figure generated successfully"
)
