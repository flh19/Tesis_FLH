#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
Genotype preprocessing, quality control, and VCF generation pipeline
===============================================================================

Description
-----------
This script performs preprocessing, quality control (QC), normalization,
and export of genotype data derived from targeted sequencing or SNP panels.

The workflow includes:

    1. Input genotype table loading
    2. Sample and variant identifier cleaning
    3. Missing value handling
    4. Detection and correction of inverted heterozygous genotypes
    5. Removal of samples with excessive missingness
    6. Removal of variants with excessive missingness
    7. Elimination of monomorphic variants
    8. Conversion of multiallelic variants into biallelic variants
    9. Hardy–Weinberg equilibrium (HWE) filtering
    10. Genotype recoding
    11. Generation of VCF output
    12. Export of cleaned genotype table

The resulting outputs are suitable for:
    - GWAS analyses,
    - Rare variant association studies,
    - Population genetics,
    - REGENIE/PLINK pipelines,
    - Downstream genomic analyses.

Input
-----
A genotype table in TSV format containing:
    - Samples as rows
    - Variants as columns
    - Genotypes encoded as REF/ALT strings

Example:
    genotype_diabet_fixed.tsv

Output
------
The script generates:
    - Filtered genotype Excel file
    - VCF file containing QC-passed variants

Main filtering steps:
    - Missingness filtering
    - Multiallelic-to-biallelic conversion
    - Hardy–Weinberg equilibrium filtering

Dependencies
------------
    pandas
    numpy
    scipy
    scikit-learn
    statsmodels
    matplotlib
    seaborn

Author
------
Prepared for genomic preprocessing and variant quality control workflows.

===============================================================================
"""

#!/usr/bin/env python3

# ============================================================
# Imports
# ============================================================

import argparse
import os
import random
import re
import warnings
from pathlib import Path
from zipfile import ZipFile

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import patsy as pt
import scipy.stats as stats
import seaborn as sns
import sklearn
import statsmodels.api as sm

from matplotlib.ticker import FuncFormatter
from scipy import stats
from scipy.stats import chi2
from scipy.stats import chisquare

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    roc_auc_score,
    roc_curve,
)

from sklearn.preprocessing import StandardScaler

from statsmodels.stats.multitest import multipletests

# Ignore warnings
warnings.filterwarnings('ignore')

# ============================================================
# Quality-control parameters
# ============================================================

# Maximum allowed missingness per sample
threshold_sample = 0.2

# Maximum allowed missingness per variant
threshold = 0.1

# ============================================================
# Load genotype data
# ============================================================

# Input genotype table
name_snp_table = (
    '/content/drive/MyDrive/share_ugdg/bioinfo/scripts/'
    'ALL/UGD_GLM_analysis/data/FISDM2_V2_SUB/'
    'genotype_diabet_fixed.tsv'
)

# Load TSV genotype table
df_snp = pd.read_csv(
    name_snp_table,
    sep='\t',
    header='infer'
)

# ============================================================
# Sample and variant identifier cleaning
# ============================================================

# Rename sample column
df_snp.rename(
    columns={'Unnamed: 0': 'sample'},
    inplace=True
)

# Use sample names as index
df_snp.index = df_snp.pop('sample')

# ------------------------------------------------------------
# Remove rs and COSMIC identifiers from column names
# ------------------------------------------------------------

for col in df_snp.columns:

    if 'rs' in col:
        new_col = col.split('-rs')[0]

    elif 'COS' in col:
        new_col = col.split('-COS')[0]

    else:
        new_col = col

    df_snp.rename(
        columns={col: new_col},
        inplace=True
    )

# ------------------------------------------------------------
# Remove A_ and _R labels from sample names
# ------------------------------------------------------------

for row in df_snp.index:

    row = str(row)

    if '_R' in row and 'A_' in row:

        new_row = (
            row
            .replace('_R', '')
            .replace('A_', '')
        )

    else:
        new_row = row

    df_snp.rename(
        index={row: new_row},
        inplace=True
    )

# ============================================================
# Missing value handling
# ============================================================

# Replace "./." with NaN
df_snp_nan = df_snp.copy()

df_snp_nan = df_snp_nan.replace(
    './.',
    np.nan
)

# ============================================================
# Detection of inverted heterozygous genotypes
# ============================================================

# Identify columns where the alternative allele appears first
columnas_alelo_alt = []

for col in df_snp_nan.columns:

    # Reference allele
    alelo_ref = col.split('-')[-2]

    # Alternative alleles
    alelos_alt = col.split('-')[-1].split(',')

    for genotipo in df_snp_nan[col].dropna():

        alelos = genotipo.split('/')

        # Ignore INDELs
        if (
            len(alelos[0]) == 1 and
            len(alelos[1]) == 1
        ):

            if (
                len(alelos) == 2 and
                alelos[0] in alelos_alt and
                alelos[1] == alelo_ref
            ):

                columnas_alelo_alt.append(col)
                break

# ============================================================
# Reorder inverted heterozygous alleles
# ============================================================

df_snp_nan_het = df_snp_nan.copy()

for col in columnas_alelo_alt:

    try:

        alelo_ref = col.split('-')[-2]
        alelos_alt = col.split('-')[-1].split(',')

    except IndexError:
        continue

    for idx, genotipo in (
        df_snp_nan_het[col]
        .dropna()
        .items()
    ):

        alelos = genotipo.split('/')

        if (
            len(alelos) == 2 and
            len(alelos[0]) == 1 and
            len(alelos[1]) == 1
        ):

            if (
                alelos[0] in alelos_alt and
                alelos[1] == alelo_ref
            ):

                # Place reference allele first
                df_snp_nan_het.at[idx, col] = (
                    f"{alelos[1]}/{alelos[0]}"
                )

# ============================================================
# Remove samples with excessive missingness
# ============================================================

missing_fraction_rows = (
    df_snp_nan_het
    .isnull()
    .mean(axis=1)
)

df_snp_nan_het_row = (
    df_snp_nan_het[
        missing_fraction_rows <= threshold_sample
    ]
    .copy()
)

print(
    f'Número inicial de muestras:',
    df_snp_nan_het.shape[0]
)

print(
    f'Número de muestras con menos de un '
    f'{threshold_sample}% de valores nulos:',
    df_snp_nan_het_row.shape[0]
)

print()

# ============================================================
# Remove variants with excessive missingness
# ============================================================

missing_fraction = (
    df_snp_nan_het_row
    .isnull()
    .mean()
)

superan_threshold = (
    missing_fraction[
        missing_fraction <= threshold
    ]
    .index
)

df_snp_nan_het_row_col = (
    df_snp_nan_het_row[
        superan_threshold
    ]
)

print(
    f'Número incial de variantes:',
    df_snp_nan_het_row.shape[1]
)

print(
    f'Número de variantes con menos de un '
    f'{threshold}% de valores nulos:',
    len(superan_threshold)
)

# ============================================================
# Remove monomorphic variants and update headers
# ============================================================

df_snp_nan_het_row_col = (
    df_snp_nan_het_row_col.copy()
)

alt_remove = []

for col in df_snp_nan_het_row_col.columns:

    print(col)

    aux = (
        df_snp_nan_het_row_col[col]
        .dropna()
        .astype(str)
        .unique()
    )

    partes = col.split('-')

    coordenadas = partes[0]
    ref_aux = partes[1]

    alt_list = []

    for combinacion in set(aux):

        alt_aux = combinacion.split('/')[1]

        if (
            alt_aux not in alt_list and
            alt_aux != ref_aux
        ):
            alt_list.append(alt_aux)

    if len(alt_list) >= 1:

        df_snp_nan_het_row_col.rename(
            columns={
                col:
                f'{coordenadas}-{ref_aux}-{",".join(alt_list)}'
            },
            inplace=True
        )

    else:
        alt_remove.append(col)

# Remove variants without alternative alleles
df_snp_nan_het_row_col_update = (
    df_snp_nan_het_row_col
    .drop(columns=alt_remove)
)

print(
    "Variantes eliminadas porque ya no contienen variantes:",
    alt_remove
)

# ============================================================
# Convert multiallelic variants into biallelic variants
# ============================================================

def genotype_frequencies(column):

    """
    Convert multiallelic variants into biallelic variants
    by retaining the most frequent alternative allele.
    """

    # Extract REF and ALT alleles
    ref = column.name.split("-")[1]

    alts = column.name.split("-")[2]

    list_alts = alts.split(',')

    # Skip already biallelic variants
    if len(list_alts) == 1:
        return column

    elif len(list_alts) > 1:

        genotypes = column.dropna()

        dict_allele = {}

        # Count allele frequencies
        for genotype in genotypes:

            alleles = genotype.split('/')

            for allele in alleles:

                dict_allele[allele] = (
                    dict_allele.get(allele, 0) + 1
                )

        # Remove REF allele
        dict_allele.pop(ref, None)

        # Keep most frequent ALT allele
        if dict_allele:

            major_allele = max(
                dict_allele,
                key=dict_allele.get
            )

            column = column.apply(
                lambda x:
                '/'.join([
                    major_allele
                    if allele in dict_allele
                    else allele
                    for allele in x.split('/')
                ])
                if pd.notna(x)
                else x
            )

    return column

# ============================================================
# Compare DataFrame columns
# ============================================================

def compare_column_contents(df1, df2):

    """
    Identify columns with different contents
    between two DataFrames.
    """

    common_columns = (
        set(df1.columns)
        .intersection(df2.columns)
    )

    different_columns = []

    for column in common_columns:

        if not df1[column].equals(df2[column]):

            different_columns.append(column)

    return different_columns

# Count multiallelic variants
cont = 0

for col in df_snp_nan_het_row_col_update.columns:

    if ',' in col:
        cont += 1

print(
    f"Number of initial multiallelic variants:{cont}"
)

# Apply multiallelic conversion
df_snp_update_filt_max_alt = (
    df_snp_nan_het_row_col_update
    .apply(genotype_frequencies)
)

# Identify modified columns
modificadas_contenido = compare_column_contents(
    df_snp_update_filt_max_alt,
    df_snp_nan_het_row_col_update
)

print("\nMultiallelic variants transformed "
      "to biallelic variants:\n")

for i in modificadas_contenido:

    print(
        "Original:\n",
        df_snp_nan_het_row_col_update[i].value_counts(),
        '\nUpdated:\n',
        df_snp_update_filt_max_alt[i].value_counts(),
        '\n'
    )

    print()

# ============================================================
# Update headers after multiallelic conversion
# ============================================================

"""
After converting multiallelic variants into biallelic variants,
the variant headers are updated again to ensure consistency
between genotype content and variant annotation.

Variants that no longer contain alternative alleles are removed.
"""

df_snp_filt_max_alt_update2 = (
    df_snp_update_filt_max_alt.copy()
)

alt_remove = []

for col in df_snp_update_filt_max_alt.columns:

    aux = (
        df_snp_update_filt_max_alt[col]
        .dropna()
        .astype(str)
        .unique()
    )

    partes = col.split('-')

    coordenadas = partes[0]
    ref_aux = partes[1]

    alt_list = []

    for combinacion in set(aux):

        alt_aux = combinacion.split('/')[1]

        if (
            alt_aux not in alt_list and
            alt_aux != ref_aux
        ):
            alt_list.append(alt_aux)

    # Update variant name
    if len(alt_list) >= 1:

        df_snp_update_filt_max_alt.rename(
            columns={
                col:
                f'{coordenadas}-{ref_aux}-{",".join(alt_list)}'
            },
            inplace=True
        )

    # Mark monomorphic variants for removal
    else:
        alt_remove.append(col)

# Remove monomorphic variants
df_snp_filt_max_alt_update2 = (
    df_snp_update_filt_max_alt
    .drop(columns=alt_remove)
)

# ============================================================
# Hardy–Weinberg equilibrium (HWE) filtering
# ============================================================

"""
Variants are tested for deviation from Hardy–Weinberg equilibrium.

Genotypes are numerically encoded:
    0 = homozygous reference
    1 = heterozygous
    2 = homozygous alternative

Variants failing Bonferroni-corrected HWE filtering
are removed from the final dataset.
"""

# ------------------------------------------------------------
# Genotype recoding
# ------------------------------------------------------------

df_numeric = df_snp_filt_max_alt_update2.copy()

for variant in df_numeric.columns:

    chr, rest = variant.split(':')

    pos, ref, alt = rest.split('-')

    dict_cod = {
        f'{ref}/{ref}': 0,
        f'{ref}/{alt}': 1,
        f'{alt}/{alt}': 2
    }

    df_numeric[variant] = (
        df_numeric[variant]
        .map(dict_cod)
    )

# ------------------------------------------------------------
# Compute HWE p-values
# ------------------------------------------------------------

hwe_p_values = []

for col in df_numeric.columns:

    genotypes = df_numeric[col].dropna()

    # Observed genotype counts
    obs_hom_ref = (genotypes == 0).sum()
    obs_het     = (genotypes == 1).sum()
    obs_hom_alt = (genotypes == 2).sum()

    total = (
        obs_hom_ref +
        obs_het +
        obs_hom_alt
    )

    # Avoid division-by-zero errors
    if total == 0:

        hwe_p_values.append(1)

        continue

    # --------------------------------------------------------
    # Allele frequencies
    # --------------------------------------------------------

    p = (
        2 * obs_hom_ref + obs_het
    ) / (2 * total)

    q = 1 - p

    # --------------------------------------------------------
    # Expected genotype counts under HWE
    # --------------------------------------------------------

    exp_hom_ref = total * p**2
    exp_het     = total * 2 * p * q
    exp_hom_alt = total * q**2

    obs = [
        obs_hom_ref,
        obs_het,
        obs_hom_alt
    ]

    exp = [
        exp_hom_ref,
        exp_het,
        exp_hom_alt
    ]

    # --------------------------------------------------------
    # Chi-square test
    # --------------------------------------------------------

    chi2_value, p_value_raw = chisquare(
        obs,
        f_exp=exp
    )

    # One degree of freedom
    p_value = 1 - chi2.cdf(
        chi2_value,
        df=1
    )

    hwe_p_values.append(p_value)

# Convert to pandas Series
hwe_p_values = pd.Series(
    hwe_p_values,
    index=df_numeric.columns
)

# ------------------------------------------------------------
# Bonferroni threshold
# ------------------------------------------------------------

hwe_threshold = (
    0.05 / df_numeric.shape[1]
)

# Keep variants passing HWE
df_hwe = df_numeric.loc[
    :,
    hwe_p_values > hwe_threshold
]

# Variants removed by HWE filtering
removed_variants = df_numeric.loc[
    :,
    hwe_p_values <= hwe_threshold
]

# Final filtered genotype table
df_snp_filt_max_alt_update2_hwe = (
    df_snp_filt_max_alt_update2[
        df_hwe.columns
    ]
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

print("Variants removed due to HWE filtering:")

print(set(removed_variants.columns))

print("HWE threshold:", hwe_threshold)

# ============================================================
# VCF generation
# ============================================================

"""
This section generates a Variant Call Format (VCF) file
from the filtered genotype matrix.

Genotypes are recoded into standard VCF notation:
    REF/REF -> 0/0
    REF/ALT -> 0/1
    ALT/ALT -> 1/1

Missing values are encoded as:
    ./.
"""

# ============================================================
# Generate VCF function
# ============================================================

def generate_vcf(df, output_file):

    """
    Generate a VCF file from genotype data.

    Parameters
    ----------
    df : pandas.DataFrame
        Genotype matrix.
    output_file : str
        Output VCF filename.
    """

    vcf_df = df.copy()

    # Sample names
    samples = vcf_df.index

    samples_str = []

    for i in samples:
        samples_str.append(str(i))

    # --------------------------------------------------------
    # VCF header
    # --------------------------------------------------------

    header_samples = (
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
        + '\t'.join(samples_str)
    )

    vcf_lines = []

    chr_list = {}

    # --------------------------------------------------------
    # Process variants
    # --------------------------------------------------------

    for variant in vcf_df.columns:

        chr = variant.split(':')[0]

        pos = (
            variant
            .split(':')[1]
            .split('-')[0]
        )

        ref = (
            variant
            .split(':')[1]
            .split('-')[1]
        )

        alt = (
            variant
            .split(':')[1]
            .split('-')[2]
        )

        # Store chromosome lengths
        if chr not in chr_list:

            chr_list[chr] = [pos]

        else:

            chr_list[chr].append(pos)

        # ----------------------------------------------------
        # Genotype recoding
        # ----------------------------------------------------

        dict_cod = {

            ref + '/' + ref: '0/0',

            ref + '/' + alt: '0/1',

            alt + '/' + alt: '1/1'
        }

        vcf_df[variant] = (
            vcf_df[variant]
            .map(dict_cod)
            .fillna("./.")
        )

        print(vcf_df)

        # Convert genotypes to string
        genotypes = '\t'.join(
            vcf_df[variant].astype(str)
        )

        # Build VCF line
        vcf_lines.append(
            f"{chr}\t{pos}\t"
            f"{str(chr)+':'+str(pos)+'-'+ref+'-'+alt}\t"
            f"{ref}\t{alt}\t.\tPASS\t.\tGT\t{genotypes}"
        )

    # ========================================================
    # Write VCF file
    # ========================================================

    with open(output_file, 'w') as f:

        # File format
        f.write('##fileformat=VCFv4.2\n')

        # Chromosome metadata
        for c in chr_list.keys():

            f.write(
                f'##contig=<ID={c},length={max(chr_list[c])}>\n'
            )

        # Source metadata
        f.write('##source=UGD\n')

        # FORMAT field
        f.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,'
            'Description="Genotype">\n'
        )

        # Sample header
        f.write(
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
            + '\t'.join(samples_str)
            + '\n'
        )

        # Variant entries
        f.write('\n'.join(vcf_lines) + '\n')

    print(f"VCF file saved to {output_file}")

# ============================================================
# Save filtered VCF file
# ============================================================

out_name_vcf = (
    name_snp_table.replace(
        '.tsv',
        '_noInvHet_filtNA_multi2bi_filtHWE_20250318.vcf'
    )
)

generate_vcf(
    df_snp_filt_max_alt_update2_hwe,
    out_name_vcf
)

# ============================================================
# Save filtered genotype table
# ============================================================

out_name = (
    name_snp_table.replace(
        '.tsv',
        '_noInvHet_filtNA_multi2bi_filtHWE_20250318.xlsx'
    )
)

df_snp_filt_max_alt_update2_hwe.to_excel(
    out_name,
    index=True
)
