#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =========================================================
# COVERAGE QC PLOTTING SCRIPT
# =========================================================
#
# Description
# -----------
# This script generates a combined Quality Control (QC)
# figure for sequencing coverage analysis.
#
# The output figure contains:
#
#   (a) A boxplot showing the distribution of coverage
#       values grouped by predefined coverage ranges.
#
#   (b) A line plot representing the percentage of bases
#       covered at different sequencing depths for all
#       samples.
#
# The script:
#   - Reads coverage statistics from an XLSX file
#   - Reads per-base coverage distributions from a CSV file
#   - Groups samples according to coverage thresholds
#   - Generates publication-quality plots
#   - Exports the final figure as PNG and PDF
#
#
# Input files
# -----------
# 1. XLSX file:
#    Must contain a numeric column with coverage values
#    (default column name: "media")
#
# 2. CSV file:
#    Expected format:
#
#       sample coverage percentage_bases
#
#    Example:
#       sample1 1 0.98
#       sample1 2 0.95
#       sample2 1 0.99
#
#
# Usage
# -----
# Example execution:
#
# python script.py \
#     --xlsx coverage_stats.xlsx \
#     --csv coverage_distribution.csv \
#     --outdir results \
#     --prefix qc_plot
#
#
# Required arguments
# ------------------
# --xlsx      Path to XLSX file containing coverage metrics
# --csv       Path to CSV file containing coverage curves
# --outdir    Output directory
#
#
# Optional arguments
# ------------------
# --prefix        Output filename prefix (default: qc)
# --value-col     Column name containing coverage values
#                  in the XLSX file (default: media)
# --ylabel        Y-axis label for the boxplot
#                  (default: Cobertura media (X))
# --threshold     Coverage threshold value
#                  (default: 20)
# --showfliers    Display boxplot outliers
#
#
# Output
# ------
# The script generates:
#
#   - <prefix>.png
#   - <prefix>.pdf
#
# inside the specified output directory.
#
# =========================================================

# =========================================================
# IMPORTS
# =========================================================
# Standard library imports
import os
import argparse

# Numerical and data processing libraries
import numpy as np
import pandas as pd

# Plotting libraries
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, MultipleLocator, FuncFormatter

# Interpolation for smooth curves
from scipy.interpolate import make_interp_spline


# =========================================================
# FORMAT UTILITIES
# =========================================================
def fmt_miles(x):
    """
    Format numbers with thin-space thousands separators.

    Example:
        10000 -> '10 000'
    """
    return f"{int(x):,}".replace(",", "\u2009")


# =========================================================
# GLOBAL PLOT STYLE
# =========================================================
def set_style():
    """
    Configure global matplotlib style parameters.

    This function standardizes:
    - Figure resolution
    - Font sizes
    - Axis appearance
    - Tick formatting
    - Background colors
    """
    plt.rcParams.update({
        "figure.dpi": 120,
        "savefig.dpi": 600,
        "font.family": "DejaVu Sans",
        "font.size": 16,
        "axes.titlesize": 16,
        "axes.labelsize": 25,
        "axes.linewidth": 1.0,
        "xtick.labelsize": 20,
        "ytick.labelsize": 20,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "grid.linewidth": 0.8,
        "axes.facecolor": "white",
        "figure.facecolor": "white",
    })


# =========================================================
# DATA LOADING FUNCTIONS
# =========================================================
def load_xlsx(xlsx_path, sheet=None):
    """
    Load XLSX file into a pandas DataFrame.

    Parameters
    ----------
    xlsx_path : str
        Path to the Excel file.
    sheet : str or int, optional
        Sheet name or index.

    Returns
    -------
    pandas.DataFrame
        Loaded dataframe.

    Raises
    ------
    ValueError
        If the XLSX file is empty.
    """
    df = pd.read_excel(xlsx_path, sheet_name=0 if sheet is None else sheet)

    if df.empty:
        raise ValueError("XLSX vacío")

    return df


def load_csv_coverage(csv_path):
    """
    Load coverage CSV file and preprocess data.

    Expected columns:
    - sample
    - coverage
    - percentage_bases

    The percentage values are converted to percentages
    multiplying by 100.

    Returns
    -------
    pandas.DataFrame
    """
    df = pd.read_csv(
        csv_path,
        sep=r"[\s,]+",
        engine="python",
        header=None,
        names=["sample", "coverage", "percentage_bases"]
    )

    # Convert columns to numeric values
    df["coverage"] = pd.to_numeric(df["coverage"], errors="coerce")
    df["percentage_bases"] = pd.to_numeric(df["percentage_bases"], errors="coerce")

    # Remove invalid rows
    df = df.dropna().copy()

    # Convert proportions to percentages
    df["percentage_bases"] *= 100

    return df


# =========================================================
# DATA VALIDATION AND GROUPING
# =========================================================
def validate_and_get_values(df, value_col):
    """
    Validate and extract numeric values from a dataframe column.

    Parameters
    ----------
    df : pandas.DataFrame
    value_col : str
        Column containing the values of interest.

    Returns
    -------
    numpy.ndarray
        Clean numeric values.
    """
    d = df.copy()

    # Convert values to numeric
    d[value_col] = pd.to_numeric(d[value_col], errors="coerce")

    # Remove NaN values
    d = d.dropna(subset=[value_col])

    return d[value_col].values


def build_groups(values, threshold):
    """
    Split coverage values into predefined coverage groups.

    Groups:
    - Total
    - ≥100X
    - 100–50X
    - 50–20X
    - <threshold X

    Returns
    -------
    list of tuples
        Each tuple contains:
        (group_name, values)
    """
    v = np.asarray(values, float)

    return [
        ("Total", v),
        ("≥100X", v[v >= 100]),
        ("100–50X", v[(v < 100) & (v >= 50)]),
        ("50–20X", v[(v < 50) & (v >= threshold)]),
        (f"<{int(threshold)}X", v[v < threshold]),
    ]


# =========================================================
# BOXPLOT FORMATTING
# =========================================================
def _format_boxplot(ax, labels, ns, ylabel, threshold):
    """
    Apply consistent formatting to the boxplot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    labels : list
        Group labels.
    ns : list
        Number of samples per group.
    ylabel : str
        Y-axis label.
    threshold : float
        Threshold value shown as a dashed horizontal line.
    """

    # Add sample counts to labels
    xticklabels = [f"{lab}\n(N={fmt_miles(n)})" for lab, n in zip(labels, ns)]

    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(xticklabels)

    ax.set_ylabel(ylabel)

    # Threshold reference line
    ax.axhline(threshold, linestyle="--", linewidth=1.6, color="#C00000")

    # Axis formatting
    ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: fmt_miles(y)))

    # Horizontal grid
    ax.grid(axis="y", alpha=0.18)


# =========================================================
# MAIN PLOTTING FUNCTION
# =========================================================
def plot_combined(groups, df_coverage, out_prefix, ylabel_box, threshold, showfliers):
    """
    Generate combined figure with:
    (a) Boxplot
    (b) Coverage lineplot

    The figure is exported in PNG and PDF formats.
    """

    # Apply global plotting style
    set_style()

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8), dpi=600)

    # Extract labels and data
    labels = [g[0] for g in groups]
    data = [g[1] for g in groups]
    ns = [len(x) for x in data]

    # =====================================================
    # BOX PLOT
    # =====================================================
    bp = ax1.boxplot(
        data,
        widths=0.52,
        patch_artist=True,
        showfliers=showfliers,

        # Median line style
        medianprops=dict(linewidth=2.2, color="#111111"),

        # Whisker style
        whiskerprops=dict(linewidth=1.3, color="#444444"),

        # Cap style
        capprops=dict(linewidth=1.3, color="#444444"),

        # Box border style
        boxprops=dict(linewidth=1.3, color="#444444"),

        # Outlier style
        flierprops=dict(
            marker="o",
            markersize=3.5,
            markerfacecolor="#444444",
            markeredgecolor="black",
            markeredgewidth=0.3,
            alpha=0.6
        ),
    )

    # Apply the same color to all boxes
    for box in bp["boxes"]:
        box.set_facecolor("#008B8B")
        box.set_alpha(0.85)

    # Apply custom formatting
    _format_boxplot(ax1, labels, ns, ylabel_box, threshold)

    # Add subplot label
    ax1.text(-0.08, 1.10, '(a)', transform=ax1.transAxes,
             fontsize=20, fontweight='bold')

    # =====================================================
    # LINE PLOT
    # =====================================================
    for sample in df_coverage['sample'].unique():

        # Select one sample and sort by coverage
        d = df_coverage[df_coverage['sample'] == sample].sort_values("coverage")

        x = d["coverage"].values
        y = d["percentage_bases"].values

        # Smooth curve if enough points are available
        if len(x) > 3:
            x_new = np.linspace(x.min(), x.max(), 300)
            y_new = make_interp_spline(x, y)(x_new)

            ax2.plot(
                x_new,
                y_new,
                alpha=0.2,
                linewidth=0.6,
                color='grey'
            )

        # Otherwise plot raw values
        else:
            ax2.plot(x, y, alpha=0.3, linewidth=0.6)

    # Axis labels (kept in Spanish intentionally)
    ax2.set_ylabel('Bases (%)')
    ax2.set_xlabel('Cobertura (X)')

    # X-axis formatting
    ax2.xaxis.set_major_formatter(FuncFormatter(lambda x, _: fmt_miles(x)))
    ax2.xaxis.set_major_locator(MultipleLocator(10))

    # Vertical threshold line
    ax2.axvline(x=20, color='red', linestyle='--', alpha=0.7)

    # Add subplot label
    ax2.text(-0.08, 1.10, '(b)', transform=ax2.transAxes,
             fontsize=20, fontweight='bold')

    # Adjust layout
    fig.tight_layout()

    # Save outputs
    fig.savefig(f"{out_prefix}.png", bbox_inches="tight")
    fig.savefig(f"{out_prefix}.pdf", bbox_inches="tight")

    # Close figure to free memory
    plt.close()


# =========================================================
# MAIN EXECUTION
# =========================================================
def main():
    """
    Main execution function.

    Workflow:
    1. Parse command-line arguments
    2. Load XLSX data
    3. Validate and group values
    4. Load coverage CSV
    5. Create output directory
    6. Generate plots
    """

    # -----------------------------------------------------
    # ARGUMENT PARSER
    # -----------------------------------------------------
    p = argparse.ArgumentParser()

    p.add_argument("--xlsx", required=True)
    p.add_argument("--csv", required=True)
    p.add_argument("--outdir", required=True)

    p.add_argument("--prefix", default="qc")
    p.add_argument("--value-col", default="media")
    p.add_argument("--ylabel", default="Cobertura media (X)")
    p.add_argument("--threshold", type=float, default=20)

    p.add_argument("--showfliers", action="store_true")

    args = p.parse_args()

    # -----------------------------------------------------
    # LOAD AND PREPARE DATA
    # -----------------------------------------------------
    df_xlsx = load_xlsx(args.xlsx)

    values = validate_and_get_values(df_xlsx, args.value_col)

    groups = build_groups(values, args.threshold)

    df_csv = load_csv_coverage(args.csv)

    # -----------------------------------------------------
    # CREATE OUTPUT DIRECTORY
    # -----------------------------------------------------
    os.makedirs(args.outdir, exist_ok=True)

    # -----------------------------------------------------
    # GENERATE PLOTS
    # -----------------------------------------------------
    plot_combined(
        groups,
        df_csv,
        os.path.join(args.outdir, args.prefix),
        args.ylabel,
        args.threshold,
        args.showfliers
    )


# =========================================================
# SCRIPT ENTRY POINT
# =========================================================
if __name__ == "__main__":
    main()
