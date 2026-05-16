#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
Principal Component Analysis (PCA) visualization and eigenvalue plots
===============================================================================

Description
-----------
This script generates publication-quality figures from principal component
analysis (PCA) results commonly obtained in population genetics and genomic
studies.

The script produces two complementary visualizations:

    1. Eigenvalue plot ("scree plot")
       - Displays the variance explained by each principal component.
       - Useful for evaluating the contribution of PCs and identifying the
         major axes of genetic variation.

    2. PCA ancestry projection
       - Displays samples projected onto PC1 and PC2.
       - Samples are colored according to ancestry classification.
       - Commonly used to assess population stratification, ancestry clustering,
         and cohort homogeneity.

These figures are typically included in:
    - GWAS studies,
    - Sequencing analyses,
    - Population genetics studies,
    - Quality control pipelines.

Input files
-----------
1. principal_components.eigenval
    File containing eigenvalues for each principal component.

2. principal_components.eigenvec
    File containing sample projections for each principal component.

3. ancestrias.txt
    Tab-delimited file containing ancestry labels for each sample.
    Required columns:
        - IID
        - ANCESTRÍA

Output
------
The script generates:
    - Eigenvalue plot:
        * varianza_explicada.png
        * varianza_explicada.pdf

    - PCA ancestry projection:
        * proyeccion_ancestria.png
        * proyeccion_ancestria.pdf

All figures are saved at 600 dpi resolution.

Usage
-----
Run directly from the command line:

    python plot_pca_ancestry.py

Dependencies
------------
    pandas
    matplotlib
    seaborn

Author
------
Prepared for population genetics and genomic cohort visualization.

===============================================================================
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter

# =========================
# General style configuration
# =========================

# Set seaborn visual style
sns.set(style="whitegrid", palette="muted")

# Function to display decimal numbers using comma separator
def comma_formatter(x, pos):
    return f'{x:.2f}'.replace('.', ',')

# =========================
# Load input data
# =========================

# Load eigenvalues
eigenval = pd.read_csv(
    'principal_components.eigenval',
    header=None,
    names=['Eigenvalue']
)

# Load eigenvectors (principal component coordinates)
eigenvec = pd.read_csv(
    'principal_components.eigenvec',
    sep='\s+',
    header=0
)

# Load ancestry labels
ancestry_data = pd.read_csv(
    'ancestrias.txt',
    sep='\t'
)

# Merge ancestry information with PCA coordinates
eigenvec = eigenvec.merge(
    ancestry_data[['IID', 'ANCESTRÍA']],
    on='IID'
)

# ============================================================
# 1. Eigenvalue plot (scree plot)
# ============================================================

# Create figure
fig, ax = plt.subplots(figsize=(10, 6), dpi=600)

# Plot eigenvalues across principal components
ax.plot(
    range(1, len(eigenval) + 1),
    eigenval['Eigenvalue'],
    marker='o',
    linestyle='-',
    color='#1f77b4',
    markersize=7
)

# Axis labels
ax.set_xlabel('Componente Principal', fontsize=12, labelpad=10)
ax.set_ylabel('Valor Propio (Eigenvalue)', fontsize=12, labelpad=10)

# Display all PC numbers on x-axis
ax.set_xticks(range(1, len(eigenval) + 1))

# Add horizontal grid
ax.grid(True, axis='y', linestyle='--', alpha=0.7)

# Use comma as decimal separator
ax.yaxis.set_major_formatter(FuncFormatter(comma_formatter))

# Customize plot borders
for spine in ax.spines.values():
    spine.set_edgecolor('black')
    spine.set_linewidth(1.5)

# Optimize layout
plt.tight_layout()

# Save figures
plt.savefig('varianza_explicada.png', dpi=600, bbox_inches='tight')
plt.savefig('varianza_explicada.pdf', dpi=600, bbox_inches='tight')

# Close figure
plt.close()

# ============================================================
# 2. PCA ancestry projection
# ============================================================

# Define color palette for ancestry groups
colores_ancestria = {
    'EAS': '#87CEFA',
    'AMR': '#8A2BE2',
    'AFR': '#E41A1C',
    'EUR': '#2E8B57',
    'SAS': '#FF8C00'
}

# Create figure
fig, ax = plt.subplots(figsize=(8.5, 4.5), dpi=600)

# Leave extra space on the right side for legend placement
fig.subplots_adjust(right=0.72)

# Generate PCA scatter plot
sns.scatterplot(
    x=eigenvec['PC1'],
    y=eigenvec['PC2'],
    hue=eigenvec['ANCESTRÍA'],
    palette=colores_ancestria,
    s=150,
    edgecolor='black',
    alpha=0.8,
    ax=ax
)

# Axis labels
ax.set_xlabel('Componente Principal 1 (PC1)', fontsize=15, labelpad=10)
ax.set_ylabel('Componente Principal 2 (PC2)', fontsize=15, labelpad=10)

# Tick label sizes
ax.tick_params(axis='both', labelsize=15)

# Add grid
ax.grid(True, linestyle='--', alpha=0.5)

# Use comma as decimal separator
ax.xaxis.set_major_formatter(FuncFormatter(comma_formatter))
ax.yaxis.set_major_formatter(FuncFormatter(comma_formatter))

# Customize plot borders
for spine in ax.spines.values():
    spine.set_edgecolor('black')
    spine.set_linewidth(1.5)

# Place legend outside the plot area on the right
ax.legend(
    title='Ancestría',
    title_fontsize=14,
    fontsize=14,
    frameon=True,
    edgecolor='black',
    fancybox=False,
    loc='center left',
    bbox_to_anchor=(1.02, 0.5)
)

# Optimize layout
plt.tight_layout()

# Save figures
plt.savefig('proyeccion_ancestria.png', dpi=600, bbox_inches='tight')
plt.savefig('proyeccion_ancestria.pdf', dpi=600, bbox_inches='tight')

# Close figure
plt.close()

# Final status message
print("Las gráficas se han guardado correctamente (PNG y PDF, DPI=600)")
