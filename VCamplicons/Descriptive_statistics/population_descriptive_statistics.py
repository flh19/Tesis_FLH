#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===========================================================================
Descriptive Statistics by T2D Status
===========================================================================

This script generates descriptive statistics stratified by Type 2 Diabetes
(T2D) status (0 = controls, 1 = cases).

OUTPUTS
-------
1. Quantitative variables:
   - Sample size
   - Mean
   - Standard deviation
   - Median
   - 25th percentile (Q1)
   - 75th percentile (Q3)
   - Interquartile range (IQR)

2. Categorical variables:
   - Total subjects in group
   - Non-missing observations
   - Counts and percentages for each category

SPECIAL HANDLING FOR SEX VARIABLE
---------------------------------
Sex coding:
   sexo = 1 → Male
   sexo = 2 → Female

FEATURES
--------
- Supports TSV and CSV files
- Supports decimal commas
- Automatically skips missing variables

=========================================================================== 
"""

# ==========================================================================
# IMPORTS
# ==========================================================================

import os
import argparse
import numpy as np
import pandas as pd


# ==========================================================================
# DEFAULT VARIABLES
# ==========================================================================

DEFAULT_QUANTITATIVE_VARS = [
    "edad",
    "IMC",
    "ICC",
    "glucosa",
    "insulina_microU_ml",
    "homa-ir",
    "CT",
    "HDL",
    "LDL",
    "TG",
    "PAS",
    "PAD"
]

DEFAULT_CATEGORICAL_VARS = [
    "HTA",
    "obesidad",
    "sexo",
    "HTC"
]


# ==========================================================================
# AUXILIARY FUNCTIONS
# ==========================================================================

def convert_to_numeric(series: pd.Series) -> pd.Series:
    """
    Convert a pandas Series to numeric format.
    Supports decimal commas.
    """

    if pd.api.types.is_numeric_dtype(series):
        return series

    s = series.astype(str).str.strip()
    s = s.str.replace(",", ".", regex=False)

    return pd.to_numeric(s, errors="coerce")


# ==========================================================================
# QUANTITATIVE VARIABLES
# ==========================================================================

def summarize_quantitative_variables(
    df,
    group_col,
    quantitative_vars
):

    results = []

    for group in sorted(df[group_col].dropna().unique()):

        df_group = df[df[group_col] == group]

        for variable in quantitative_vars:

            if variable not in df_group.columns:
                continue

            values = convert_to_numeric(
                df_group[variable]
            ).dropna().astype(float)

            n = int(values.shape[0])

            if n == 0:

                results.append({
                    "Variable": variable,
                    "T2D_Status": int(group),
                    "N_Individuals": 0,
                    "Mean": np.nan,
                    "Standard_Deviation": np.nan,
                    "Median": np.nan,
                    "Percentile_25": np.nan,
                    "Percentile_75": np.nan,
                    "Interquartile_Range": np.nan
                })

                continue

            q1 = float(values.quantile(0.25))
            q3 = float(values.quantile(0.75))

            results.append({
                "Variable": variable,
                "T2D_Status": int(group),
                "N_Individuals": n,
                "Mean": float(values.mean()),
                "Standard_Deviation": (
                    float(values.std(ddof=1))
                    if n > 1 else np.nan
                ),
                "Median": float(values.median()),
                "Percentile_25": q1,
                "Percentile_75": q3,
                "Interquartile_Range": float(q3 - q1)
            })

    return (
        pd.DataFrame(results)
        .sort_values(["Variable", "T2D_Status"])
    )


# ==========================================================================
# CATEGORICAL VARIABLES
# ==========================================================================

def summarize_categorical_variables(
    df,
    group_col,
    categorical_vars
):

    results = []

    for group in sorted(df[group_col].dropna().unique()):

        df_group = df[df[group_col] == group]

        total_group_size = int(df_group.shape[0])

        for variable in categorical_vars:

            if variable not in df_group.columns:
                continue

            values = pd.to_numeric(
                df_group[variable],
                errors="coerce"
            )

            non_missing = int(values.notna().sum())

            # ==============================================================
            # SEX VARIABLE
            # ==============================================================

            if variable.lower() == "sexo":

                # 1 = Male
                # 2 = Female

                count_1 = int((values == 1).sum())
                count_2 = int((values == 2).sum())

                pct_1 = (
                    count_1 / non_missing * 100
                    if non_missing > 0 else np.nan
                )

                pct_2 = (
                    count_2 / non_missing * 100
                    if non_missing > 0 else np.nan
                )

                category_1_label = "Male"
                category_2_label = "Female"

            # ==============================================================
            # STANDARD BINARY VARIABLES
            # ==============================================================

            else:

                count_1 = int((values == 1).sum())
                count_2 = int((values == 0).sum())

                pct_1 = (
                    count_1 / non_missing * 100
                    if non_missing > 0 else np.nan
                )

                pct_2 = (
                    count_2 / non_missing * 100
                    if non_missing > 0 else np.nan
                )

                category_1_label = "Value_1"
                category_2_label = "Value_0"

            results.append({

                "Variable": variable,
                "T2D_Status": int(group),

                "Total_Group_Size": total_group_size,
                "Non_Missing_Observations": non_missing,

                "Category_1_Label": category_1_label,
                "Category_1_Count": count_1,
                "Category_1_Percentage": pct_1,

                "Category_2_Label": category_2_label,
                "Category_2_Count": count_2,
                "Category_2_Percentage": pct_2
            })

    return (
        pd.DataFrame(results)
        .sort_values(["Variable", "T2D_Status"])
    )


# ==========================================================================
# MAIN
# ==========================================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--infile",
        required=True
    )

    parser.add_argument(
        "--sep",
        default="\t"
    )

    parser.add_argument(
        "--outdir",
        default="descriptive_statistics_output"
    )

    parser.add_argument(
        "--group-col",
        default="DM2"
    )

    parser.add_argument(
        "--quant",
        nargs="*",
        default=DEFAULT_QUANTITATIVE_VARS
    )

    parser.add_argument(
        "--cat",
        nargs="*",
        default=DEFAULT_CATEGORICAL_VARS
    )

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ----------------------------------------------------------------------
    # LOAD DATA
    # ----------------------------------------------------------------------

    df = pd.read_csv(
        args.infile,
        sep=args.sep,
        dtype=str
    )

    if args.group_col not in df.columns:

        raise ValueError(
            f"Grouping column '{args.group_col}' not found."
        )

    df[args.group_col] = pd.to_numeric(
        df[args.group_col],
        errors="coerce"
    )

    # ----------------------------------------------------------------------
    # QUANTITATIVE VARIABLES
    # ----------------------------------------------------------------------

    quantitative_results = summarize_quantitative_variables(
        df,
        args.group_col,
        args.quant
    )

    quantitative_output = os.path.join(
        args.outdir,
        "quantitative_descriptives_by_T2D.csv"
    )

    quantitative_results.to_csv(
        quantitative_output,
        index=False
    )

    # ----------------------------------------------------------------------
    # CATEGORICAL VARIABLES
    # ----------------------------------------------------------------------

    categorical_results = summarize_categorical_variables(
        df,
        args.group_col,
        args.cat
    )

    categorical_output = os.path.join(
        args.outdir,
        "categorical_descriptives_by_T2D.csv"
    )

    categorical_results.to_csv(
        categorical_output,
        index=False
    )

    # ----------------------------------------------------------------------
    # REPORT
    # ----------------------------------------------------------------------

    print("\n====================================================")
    print("Descriptive statistics successfully generated")
    print("====================================================")

    print(f"\nQuantitative output:")
    print(f"  {quantitative_output}")

    print(f"\nCategorical output:")
    print(f"  {categorical_output}")


# ==========================================================================
# SCRIPT EXECUTION
# ==========================================================================

if __name__ == "__main__":
    main()
