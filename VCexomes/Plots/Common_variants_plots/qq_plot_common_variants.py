"""
============================================================
Combined QQ Plot Generator for GWAS Association Results
============================================================

This script generates a combined Quantile-Quantile (QQ) plot
for multiple phenotypes using REGENIE association results.

"""

# ============================================================
# IMPORTS
# ============================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import chi2
from matplotlib.ticker import FuncFormatter

import os


# ============================================================
# GLOBAL STYLE CONFIGURATION
# ============================================================

plt.rcParams.update({

    "font.family": "sans-serif",

    "font.size": 14,

    "axes.linewidth": 1.1,
    "axes.labelsize": 16,

    "xtick.labelsize": 14,
    "ytick.labelsize": 14,

    "legend.fontsize": 11,

    "pdf.fonttype": 42
})


# ============================================================
# INPUT FILES
# ============================================================

# Mapping between phenotype names and REGENIE result files

FILE_PATHS = {

    "DM2": "comunes_allpopulation_DM2_DM2.regenie",

    "Glucosa": "comunes_allpopulation_glucosa_glucosa.regenie",

    "Insulina": "comunes_allpopulation_insulina_insulina.regenie",

    "HOMA-IR": "comunes_allpopulation_homa-ir_homa-ir.regenie",
}


# ============================================================
# PHENOTYPE COLORS
# ============================================================

COLORS = {

    "DM2": "#8A7CA8",

    "Glucosa": "#5B7C99",

    "Insulina": "#7A9E7E",

    "HOMA-IR": "#B07D62",
}


# ============================================================
# AXIS LABEL FORMATTER
# ============================================================

def comma_formatter(x, pos):
    """
    Format decimal numbers using commas instead of dots.

    Example:
        2.5 -> 2,5
    """

    if abs(x - int(x)) < 1e-9:
        return f"{int(x)}"

    return f"{x:.1f}".replace('.', ',')


# ============================================================
# GENOMIC INFLATION FACTOR (LAMBDA GC)
# ============================================================

def calculate_lambda(p_values):
    """
    Compute the genomic inflation factor (Lambda GC).

    Lambda GC evaluates systematic inflation in GWAS test
    statistics, often caused by population stratification
    or cryptic relatedness.

    Parameters
    ----------
    p_values : array-like
        Array of p-values.

    Returns
    -------
    float
        Genomic inflation factor.
    """

    p_values = np.array(p_values)

    # --------------------------------------------------------
    # Remove NaN values
    # --------------------------------------------------------

    p_values = p_values[~np.isnan(p_values)]

    # --------------------------------------------------------
    # Remove invalid p-values
    # --------------------------------------------------------

    p_values = p_values[(p_values > 0) & (p_values < 1)]

    if len(p_values) == 0:
        return np.nan

    # Convert p-values to chi-square statistics
    chi2_stats = chi2.ppf(1 - p_values, df=1)

    # Lambda GC
    return np.nanmedian(chi2_stats) / chi2.ppf(0.5, df=1)


# ============================================================
# LOAD AND CLEAN P-VALUES
# ============================================================

def load_clean_data_raw_p(path):
    """
    Load and clean p-values from a REGENIE result file.

    Parameters
    ----------
    path : str
        Path to GWAS results file.

    Returns
    -------
    numpy.ndarray or None
        Array of valid p-values or None if loading fails.
    """

    # --------------------------------------------------------
    # Check file existence
    # --------------------------------------------------------

    if not os.path.exists(path):

        print(f"File not found: {path}")

        return None

    try:

        # Read whitespace-separated file
        df = pd.read_csv(path, sep=r'\s+')

        # Ensure p-value column exists
        if 'Pval' not in df.columns:

            print(f"Pval column not found in {path}")

            return None

        # Convert to numeric values
        p_vals = pd.to_numeric(
            df['Pval'],
            errors='coerce'
        ).dropna().values

        # Keep only valid p-values
        return p_vals[(p_vals > 0) & (p_vals <= 1)]

    except Exception as e:

        print(f"Error reading {path}: {e}")

        return None


# ============================================================
# QQ PLOT GENERATION
# ============================================================

def draw_qq_all(ax):
    """
    Draw a combined QQ plot for all phenotypes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axis.
    """

    max_obs = 0
    max_exp = 0

    # --------------------------------------------------------
    # Null distribution reference line
    # --------------------------------------------------------

    null_line, = ax.plot(

        [0, 20],
        [0, 20],

        color='black',

        linestyle='--',

        lw=1.2,

        zorder=1
    )

    # --------------------------------------------------------
    # Plot each phenotype
    # --------------------------------------------------------

    for pheno in FILE_PATHS.keys():

        p_raw = load_clean_data_raw_p(FILE_PATHS[pheno])

        if p_raw is None:
            continue

        # Compute Lambda GC
        lam = calculate_lambda(p_raw)

        # Observed distribution
        obs_sorted = np.sort(
            -np.log10(p_raw)
        )[::-1]

        n = len(obs_sorted)

        # Expected null distribution
        exp = -np.log10(
            np.arange(1, n + 1) / (n + 1)
        )

        # Scatter plot
        ax.scatter(

            exp,
            obs_sorted,

            s=20,

            color=COLORS[pheno],

            label=(
                f"{pheno} "
                f"($\\lambda = "
                f"{str(f'{lam:.2f}').replace('.', ',')}$)"
            ),

            alpha=0.6,

            edgecolors='none',

            zorder=2,

            rasterized=True
        )

        # Update axis limits
        max_obs = max(max_obs, np.max(obs_sorted))
        max_exp = max(max_exp, np.max(exp))

    # --------------------------------------------------------
    # Axis labels
    # --------------------------------------------------------

    ax.set_xlabel(r"$-\log_{10}(p)$ esperado")

    ax.set_ylabel(r"$-\log_{10}(p)$ observado")

    # --------------------------------------------------------
    # Decimal formatting
    # --------------------------------------------------------

    ax.xaxis.set_major_formatter(
        FuncFormatter(comma_formatter)
    )

    ax.yaxis.set_major_formatter(
        FuncFormatter(comma_formatter)
    )

    # --------------------------------------------------------
    # Main legend
    # --------------------------------------------------------

    ax.legend(
        loc='upper left',
        frameon=False,
        fontsize=11
    )

    # --------------------------------------------------------
    # Reference line legend
    # --------------------------------------------------------

    from matplotlib.legend import Legend

    leg_ref = Legend(

        ax,

        [null_line],

        ['Distribución nula'],

        loc='lower right',

        frameon=False,

        fontsize=11,

        handlelength=2
    )

    ax.add_artist(leg_ref)

    # --------------------------------------------------------
    # Background grid
    # --------------------------------------------------------

    ax.grid(
        True,
        linestyle='--',
        alpha=0.3,
        zorder=0
    )

    # --------------------------------------------------------
    # Axis limits
    # --------------------------------------------------------

    ax.set_xlim(0, max_exp + 0.2)

    ax.set_ylim(0, max_obs + 1.0)


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    """
    Main execution function.
    """

    # --------------------------------------------------------
    # Create figure
    # --------------------------------------------------------

    fig, ax = plt.subplots(figsize=(5.2, 5.2))

    # Generate QQ plot
    draw_qq_all(ax)

    plt.tight_layout()

    # --------------------------------------------------------
    # Save PNG
    # --------------------------------------------------------

    plt.savefig(
        "qq_plot_common_variants.png",
        dpi=400,
        bbox_inches='tight',
        pad_inches=0.02
    )

    # --------------------------------------------------------
    # Save PDF
    # --------------------------------------------------------

    plt.savefig(
        "qq_plot_common_variants.pdf",
        dpi=300,
        bbox_inches='tight',
        pad_inches=0.02
    )

    print("Combined QQ plot generated successfully.")


# ============================================================
# SCRIPT ENTRY POINT
# ============================================================

if __name__ == "__main__":

    main()
