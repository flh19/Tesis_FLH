# Fully Commented Genotype–Phenotype Barplot Generator

```python
#!/usr/bin/env python3

"""
============================================================
GENOTYPE-PHENOTYPE BARPLOT GENERATOR
============================================================

Description:
    This script generates publication-quality barplots showing
    the distribution of quantitative phenotypes across genotype
    groups for selected genetic variants.

    The script supports:
        - Individual variant plots
        - Multipanel figures
        - Automatic genotype label generation
        - Standard deviation error bars
        - Optional exclusion of treated individuals
        - Automatic formatting of phenotype labels

    Figures are exported in both PNG and PDF formats.

------------------------------------------------------------
INPUT FILE REQUIREMENTS
------------------------------------------------------------

1. DATA FILE
------------------------------------------------------------

The main input table must contain:
    - Sample identifier column
    - Genotype columns
    - Phenotype columns
    - Covariate columns

Accepted sample ID column names:
    - IID
    - FID
    - sample_id
    - ID

Example:

    IID     rs12345     HOMA_IR     age     sex
    S1      0           1.45        54      1
    S2      1           2.12        60      0

------------------------------------------------------------
2. CONFIGURATION FILE
------------------------------------------------------------

Each line must contain:

    variant phenotype covariates filter_file model rsid gene

Example:

    chr1:12345:A:G HOMA_IR age,sex treatment.txt ADD rs12345 GCKR

Fields:
    1. Variant column name
    2. Phenotype column
    3. Covariates (comma-separated)
    4-5. Optional filter files or model information
    6. rsID label for plotting
    7. Gene symbol

------------------------------------------------------------
USAGE
------------------------------------------------------------

Basic execution:

    python script.py

Or modify the final function call:

    plot_variants_barplot(
        'datos_diabet.txt',
        'configuracion.txt'
    )

------------------------------------------------------------
OUTPUT
------------------------------------------------------------

For each variant:

    - barplot_<variant>_<phenotype>.png
    - barplot_<variant>_<phenotype>.pdf

Additionally:

    - barplot_multipanel.png
    - barplot_multipanel.pdf

------------------------------------------------------------
DEPENDENCIES
------------------------------------------------------------

Required Python packages:

    - pandas
    - numpy
    - matplotlib
    - seaborn

Install with:

    pip install pandas numpy matplotlib seaborn

============================================================
"""

# ============================================================
# IMPORTS
# ============================================================

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import numpy as np

from matplotlib.backends.backend_pdf import PdfPages


# ============================================================
# GENOTYPE LABEL GENERATION
# ============================================================


def generate_genotype_labels(variant_id):
    """
    Generate genotype labels from a variant identifier.

    Example:
        chr1-12345-A-G

    Returns:
        {
            0: "AA",
            1: "AG",
            2: "GG"
        }
    """

    try:
        parts = variant_id.split('-')

        ref, alt = parts[-2], parts[-1]

        return {
            0: f"{ref}{ref}",
            1: f"{ref}{alt}",
            2: f"{alt}{alt}"
        }

    except:
        return {
            0: "0",
            1: "1",
            2: "2"
        }


# ============================================================
# PHENOTYPE LABEL FORMATTING
# ============================================================


def format_phenotype_label(phenotype):
    """
    Convert phenotype names into cleaner plot labels.
    """

    phenotype_lower = phenotype.lower()

    if phenotype_lower == 'homa_ir' or phenotype_lower == 'homa-ir':
        return 'HOMA-IR'

    if phenotype_lower == 'log_tg':
        return 'Log(TG)'

    return phenotype.replace('_', ' ').title()


# ============================================================
# THOUSAND SEPARATOR FORMATTER
# ============================================================


def format_with_thin_space(x, pos):
    """
    Format numeric values using Unicode thin spaces as
    thousand separators for improved figure aesthetics.
    """

    thin_space = '\u2009'

    if float(x).is_integer():

        return f"{int(x):,}".replace(",", thin_space)

    else:

        integer_part, decimal_part = (
            f"{x:g}".split('.')
            if '.' in f"{x:g}"
            else (f"{x:g}", None)
        )

        integer_part = (
            f"{int(integer_part):,}"
            .replace(",", thin_space)
        )

        if decimal_part:
            return f"{integer_part},{decimal_part}"
        else:
            return integer_part


# ============================================================
# SINGLE BARPLOT CREATION
# ============================================================


def create_single_barplot(
    ax,
    subset,
    variant,
    phenotype,
    rsid,
    order,
    labels_map,
    colors,
    panel_label=None,
    gene=None
):
    """
    Generate a single genotype-phenotype barplot.

    The plot includes:
        - Mean phenotype values
        - Standard deviation error bars
        - Genotype sample sizes
        - Optional panel labels
        - Optional gene names
    """

    stats_data = []

    new_xticklabels = []

    # --------------------------------------------------------
    # SUMMARY STATISTICS CALCULATION
    # --------------------------------------------------------

    for label in order:

        group_data = subset[
            subset['Genotype_Label'] == label
        ][phenotype]

        n_size = len(group_data)

        mean_val = group_data.mean()

        std_val = group_data.std()

        stats_data.append({
            'Genotype_Label': label,
            'mean': mean_val,
            'std': std_val,
            'n': n_size
        })

        new_xticklabels.append(
            f"{label}\n(N = {format_with_thin_space(n_size, None)})"
        )

    stats_df = pd.DataFrame(stats_data)

    x_positions = np.arange(len(order))

    # --------------------------------------------------------
    # BARPLOT
    # --------------------------------------------------------

    bars = ax.bar(
        x_positions,
        stats_df['mean'],
        width=0.6,
        color=colors,
        edgecolor='black',
        linewidth=1.5,
        alpha=0.75
    )

    # --------------------------------------------------------
    # ERROR BARS
    # --------------------------------------------------------

    ax.errorbar(
        x_positions,
        stats_df['mean'],
        yerr=stats_df['std'],
        fmt='none',
        ecolor='black',
        elinewidth=1.5,
        capsize=5,
        capthick=1.5,
        zorder=3
    )

    # --------------------------------------------------------
    # X-AXIS FORMATTING
    # --------------------------------------------------------

    ax.set_xticks(x_positions)

    ax.set_xticklabels(
        new_xticklabels,
        fontsize=20
    )

    # --------------------------------------------------------
    # Y-AXIS FORMATTING
    # --------------------------------------------------------

    ax.tick_params(axis='y', labelsize=20)

    ax.yaxis.set_major_formatter(
        ticker.FuncFormatter(format_with_thin_space)
    )

    # --------------------------------------------------------
    # PANEL LABEL
    # --------------------------------------------------------

    if panel_label:

        ax.text(
            -0.15,
            1.05,
            f"({panel_label})",
            transform=ax.transAxes,
            fontsize=20,
            fontweight='bold',
            va='top',
            ha='left'
        )

    # --------------------------------------------------------
    # X-AXIS LABEL
    # --------------------------------------------------------

    if gene:

        ax.set_xlabel(
            rf"{rsid} ($\it{{{gene}}}$)",
            fontsize=20,
            labelpad=12
        )

    else:

        ax.set_xlabel(
            f"{rsid}",
            fontsize=20
        )

    # --------------------------------------------------------
    # Y-AXIS LABEL
    # --------------------------------------------------------

    ax.set_ylabel(
        format_phenotype_label(phenotype),
        fontsize=20
    )

    # --------------------------------------------------------
    # Y-LIMITS
    # --------------------------------------------------------

    y_max = ax.get_ylim()[1]

    ax.set_ylim(0, y_max)

    # --------------------------------------------------------
    # AESTHETICS
    # --------------------------------------------------------

    sns.despine(ax=ax, trim=True)

    ax.yaxis.grid(
        True,
        linestyle='--',
        alpha=0.3,
        linewidth=0.5
    )

    ax.set_axisbelow(True)


# ============================================================
# MAIN PLOTTING FUNCTION
# ============================================================


def plot_variants_barplot(data_path, config_path):
    """
    Generate all individual and multipanel barplots
    defined in the configuration file.
    """

    # --------------------------------------------------------
    # LOAD MAIN DATASET
    # --------------------------------------------------------

    df = pd.read_csv(
        data_path,
        sep=None,
        engine='python'
    )

    # --------------------------------------------------------
    # DETECT SAMPLE ID COLUMN
    # --------------------------------------------------------

    id_col = next(
        (
            c for c in [
                'IID',
                'FID',
                'sample_id',
                'ID'
            ]
            if c in df.columns
        ),
        df.columns[0]
    )

    # --------------------------------------------------------
    # LOAD CONFIGURATION FILE
    # --------------------------------------------------------

    with open(config_path, 'r') as f:

        lines = [
            line.strip()
            for line in f
            if line.strip()
        ]

    # --------------------------------------------------------
    # GLOBAL PLOT SETTINGS
    # --------------------------------------------------------

    sns.set_style("ticks")

    plt.rcParams['font.family'] = 'Arial'

    plt.rcParams['font.size'] = 14

    # Color palette for genotype groups
    colors = [
        '#4682B4',
        '#B0BEC5',
        '#7A8F7B'
    ]

    plot_data_list = []

    # ========================================================
    # PROCESS EACH CONFIGURATION ENTRY
    # ========================================================

    for line in lines:

        parts = line.split()

        variant = parts[0]

        phenotype = parts[1]

        covariates = parts[2].split(',')

        rsid = parts[5] if len(parts) >= 6 else variant

        gene = parts[6] if len(parts) >= 7 else None

        filter_file = next(
            (
                p for p in parts[3:5]
                if p.endswith('.txt')
            ),
            None
        )

        print(f"\nProcesando {rsid} ({variant})...")

        # ----------------------------------------------------
        # SUBSET REQUIRED COLUMNS
        # ----------------------------------------------------

        cols_needed = [
            id_col,
            variant,
            phenotype
        ] + covariates

        subset = df[cols_needed].copy()

        # ----------------------------------------------------
        # OPTIONAL SAMPLE EXCLUSION
        # ----------------------------------------------------

        if filter_file and os.path.exists(filter_file):

            tto_df = pd.read_csv(
                filter_file,
                sep=None,
                engine='python'
            )

            ids_a_excluir = (
                tto_df['IID']
                .astype(str)
                .unique()
            )

            subset = subset[
                ~subset[id_col]
                .astype(str)
                .isin(ids_a_excluir)
            ]

        # ----------------------------------------------------
        # REMOVE MISSING VALUES
        # ----------------------------------------------------

        subset = subset.dropna(
            subset=[variant, phenotype] + covariates
        )

        if subset.empty:

            print(
                f"  [!] Sin datos suficientes para {rsid}"
            )

            continue

        # ----------------------------------------------------
        # GENOTYPE LABELS
        # ----------------------------------------------------

        labels_map = generate_genotype_labels(variant)

        subset['Genotype_Label'] = (
            subset[variant]
            .astype(int)
            .map(labels_map)
        )

        order = [
            labels_map[0],
            labels_map[1],
            labels_map[2]
        ]

        # ----------------------------------------------------
        # STORE DATA FOR MULTIPANEL FIGURE
        # ----------------------------------------------------

        plot_data_list.append({
            'subset': subset,
            'variant': variant,
            'phenotype': phenotype,
            'rsid': rsid,
            'order': order,
            'labels_map': labels_map,
            'gene': gene
        })

        # ====================================================
        # INDIVIDUAL PLOT
        # ====================================================

        fig, ax = plt.subplots(
            figsize=(8, 6)
        )

        create_single_barplot(
            ax,
            subset,
            variant,
            phenotype,
            rsid,
            order,
            labels_map,
            colors,
            gene=gene
        )

        # ----------------------------------------------------
        # OUTPUT FILE NAME
        # ----------------------------------------------------

        clean_name = (
            rsid
            if rsid.startswith("rs")
            else variant.replace(':', '_')
        )

        plt.tight_layout()

        # ----------------------------------------------------
        # SAVE PNG
        # ----------------------------------------------------

        plt.savefig(
            f"barplot_{clean_name}_{phenotype}.png",
            dpi=600,
            bbox_inches='tight',
            facecolor='white'
        )

        print(
            f"  [OK] PNG guardado: "
            f"barplot_{clean_name}_{phenotype}.png"
        )

        # ----------------------------------------------------
        # SAVE PDF
        # ----------------------------------------------------

        plt.savefig(
            f"barplot_{clean_name}_{phenotype}.pdf",
            dpi=600,
            bbox_inches='tight',
            facecolor='white'
        )

        print(
            f"  [OK] PDF guardado: "
            f"barplot_{clean_name}_{phenotype}.pdf"
        )

        plt.close()

    # ========================================================
    # MULTIPANEL FIGURE
    # ========================================================

    if len(plot_data_list) > 0:

        print(
            f"\n--- Creando multipanel con "
            f"{len(plot_data_list)} gráficos ---"
        )

        n_plots = len(plot_data_list)

        n_cols = 2

        n_rows = int(np.ceil(n_plots / n_cols))

        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(16, 7 * n_rows)
        )

        if n_plots == 1:
            axes = np.array([axes])

        axes = axes.flatten()

        # Panel labels: a, b, c, ...
        panel_labels = [
            chr(97 + i)
            for i in range(n_plots)
        ]

        # ----------------------------------------------------
        # GENERATE EACH PANEL
        # ----------------------------------------------------

        for idx, plot_data in enumerate(plot_data_list):

            create_single_barplot(
                axes[idx],
                plot_data['subset'],
                plot_data['variant'],
                plot_data['phenotype'],
                plot_data['rsid'],
                plot_data['order'],
                plot_data['labels_map'],
                colors,
                panel_label=panel_labels[idx],
                gene=plot_data['gene']
            )

        # ----------------------------------------------------
        # HIDE UNUSED PANELS
        # ----------------------------------------------------

        for idx in range(n_plots, len(axes)):
            axes[idx].axis('off')

        # ----------------------------------------------------
        # PANEL SPACING
        # ----------------------------------------------------

        plt.subplots_adjust(
            hspace=0.4,
            wspace=0.3
        )

        # ----------------------------------------------------
        # SAVE MULTIPANEL PNG
        # ----------------------------------------------------

        plt.savefig(
            "barplot_multipanel.png",
            dpi=600,
            bbox_inches='tight',
            facecolor='white'
        )

        print(
            f"  [OK] Multipanel PNG guardado: "
            f"barplot_multipanel.png"
        )

        # ----------------------------------------------------
        # SAVE MULTIPANEL PDF
        # ----------------------------------------------------

        plt.savefig(
            "barplot_multipanel.pdf",
            dpi=600,
            bbox_inches='tight',
            facecolor='white'
        )

        print(
            f"  [OK] Multipanel PDF guardado: "
            f"barplot_multipanel.pdf"
        )

        plt.close()

    # ========================================================
    # FINAL MESSAGE
    # ========================================================

    print("\n¡Proceso completado!")


# ============================================================
# SCRIPT EXECUTION
# ============================================================

plot_variants_barplot(
    'datos_diabet.txt',
    'configuracion.txt'
)
```
