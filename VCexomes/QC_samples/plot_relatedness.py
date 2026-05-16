#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
Relatedness coefficient (φ) distribution plot
===============================================================================

Description
-----------
This script generates a publication-quality histogram showing the distribution
of pairwise relatedness coefficients (RELATEDNESS_PHI) across all sample pairs
from a vcftools/relatedness2 output file.

The figure includes:
    1. A histogram of pairwise φ values between different individuals
       (non self-pairs).
    2. A highlighted bar representing self-identity pairs.
    3. A vertical dashed line indicating a user-defined relatedness threshold.

The y-axis is displayed in logarithmic scale to improve visualization of the
distribution across multiple orders of magnitude.

This type of figure is commonly used in genetic studies and sequencing cohorts
to:
    - Identify related individuals,
    - Visualize cryptic relatedness,
    - Define sample exclusion thresholds,
    - Assess sample identity quality.

Input
-----
A tab-delimited file containing at least the following columns:
    - INDV1
    - INDV2
    - RELATEDNESS_PHI

The file is expected to be generated from vcftools relatedness analysis.

Output
------
The script generates:
    - PNG figure (600 dpi)
    - PDF figure

Usage
-----
Run directly from the command line:

    python plot_relatedness_phi.py

Dependencies
------------
    pandas
    matplotlib
    numpy

Author
------
Prepared for genetic association / sequencing cohort visualization.

===============================================================================
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter

# =========================
# 1. Load data
# =========================

df = pd.read_csv(
    "VCexomes_raw_allSamples_snps_indels_ann_sort_snv_biallelic_filter_tags_ID_HWE1e8-vcftools_relatedness2.relatedness2-edit",
    sep="\t"
)

# =========================
# 2. Separate self-pairs and non self-pairs
# =========================

# Self-pairs correspond to comparisons of each sample with itself
df_self = df[df["INDV1"] == df["INDV2"]].copy()

# Non self-pairs correspond to comparisons between different individuals
df_no_self = df[df["INDV1"] != df["INDV2"]].copy()

# =========================
# 3. Parameters
# =========================

# Relatedness threshold (commonly used for duplicate/related sample filtering)
umbral_phi = 0.33

# Total number of unique samples in the cohort
n_muestras_total = df_self["INDV1"].nunique()

# =========================
# 4. Histogram of pairwise relatedness
# =========================

# Remove missing values from phi distribution
phi_no_self = df_no_self["RELATEDNESS_PHI"].dropna()

# Define histogram bins
bins = np.linspace(phi_no_self.min(), phi_no_self.max(), 80)

# Create figure
fig, ax = plt.subplots(figsize=(7.5, 3.8), dpi=600)

# Leave extra space on the right side for the legend
fig.subplots_adjust(right=0.75)

# Plot histogram for non self-pairs
ax.hist(
    phi_no_self,
    bins=bins,
    color="lightgrey",
    edgecolor="black",
    linewidth=0.3,
    label="Pares entre muestras",
    zorder=1
)

# Use logarithmic scale for y-axis
ax.set_yscale("log")

# =========================
# 5. Self-identity bar
# =========================

# Position of self-identity bar
x_self = 0.5

# Current upper y-axis limit
y_max = ax.get_ylim()[1]

# Unicode thin space separator (U+2009) for thousands formatting
thin_space = "\u2009"

# Format total sample count
n_muestras_total_fmt = f"{n_muestras_total:,}".replace(",", thin_space)

# Plot self-identity bar
ax.bar(
    x_self,
    y_max * 0.9,
    width=0.01,
    color="#008B8B",
    edgecolor="black",
    linewidth=0.5,
    zorder=3,
    label=f"Self-identity (N = {n_muestras_total_fmt})"
)

# =========================
# 6. φ threshold line
# =========================

# Add vertical dashed threshold line
ax.axvline(
    umbral_phi,
    color="darkred",
    linestyle="--",
    linewidth=1,
    zorder=2
)

# Add threshold label
ax.text(
    umbral_phi + 0.005,
    y_max * 0.95,
    "φ = 0,33",
    fontsize=10,
    ha="left",
    va="top"
)

# =========================
# 7. Figure aesthetics
# =========================

# Axis labels
ax.set_xlabel("Coeficiente de parentesco (φ)", fontsize=11)
ax.set_ylabel("Número de pares (log)", fontsize=11)

# Tick label sizes
ax.tick_params(axis="x", labelsize=11)
ax.tick_params(axis="y", labelsize=11)

# X-axis limits
ax.set_xlim(phi_no_self.min(), 0.52)

# Remove upper and right borders
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# =========================
# 8. Legend
# =========================

# Place legend outside the plot area on the right
ax.legend(
    frameon=False,
    fontsize=10,
    loc="center left",
    bbox_to_anchor=(1.01, 0.05)
)

# =========================
# 9. X-axis formatting
# =========================

# Use comma as decimal separator
ax.xaxis.set_major_formatter(
    FuncFormatter(lambda x, _: str(round(x, 2)).replace('.', ','))
)

# Optimize layout
plt.tight_layout()

# =========================
# 10. Save figure
# =========================

# Save high-resolution PNG
plt.savefig(
    "distribucion_coeficiente_parentesco_phi_final.png",
    dpi=600,
    bbox_inches="tight"
)

# Save vector PDF
plt.savefig(
    "distribucion_coeficiente_parentesco_phi_final.pdf",
    bbox_inches="tight"
)

# Close figure
plt.close()
