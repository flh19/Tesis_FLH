#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===========================================================================
Statistical Tests for Baseline Characteristics by T2D Status
===========================================================================

This script calculates p-values and statistical tests for comparisons
between controls (T2D = 0) and cases (T2D = 1).

The statistical test applied depends on the predefined variable type.

STATISTICAL TESTS
-----------------

1. Quantitative variables reported as Mean ± SD:
   → Welch's t-test
     (Student t-test with unequal variances)

2. Quantitative variables reported as Median [IQR]:
   → Mann–Whitney U test (two-sided)

3. Categorical variables:
   → Chi-square test
   → Fisher's exact test automatically applied
     when expected frequencies < 5 in 2x2 tables

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

OUTPUT FILES
------------
1. statistical_tests_results.csv
   Contains:
   - Variable
   - Variable type
   - Statistical test used
   - P-value
   - Sample size per group

2. contingency_tables.csv
   Contains contingency tables for categorical variables
   for traceability and reproducibility.

USAGE
-----
python3 statistical_tests_by_T2D.py \
    --infile data.tsv \
    --outdir results \
    --sep "\t"

=========================================================================== 
"""

# ==========================================================================
# IMPORTS
# ==========================================================================

import os
import argparse
import numpy as np
import pandas as pd
from scipy import stats


# ==========================================================================
# USER-DEFINED VARIABLE CATEGORIES
# ==========================================================================

# Variables summarized as Mean ± SD
MEAN_SD_VARIABLES = [
    "edad",
    "IMC",
    "CT",
    "HDL",
    "LDL",
    "PAS",
    "PAD",
    "ICC"
]

# Variables summarized as Median [IQR]
MEDIAN_IQR_VARIABLES = [
    "glucosa",
    "insulina_microU_ml",
    "homa-ir",
    "TG"
]

# Binary or categorical variables
CATEGORICAL_VARIABLES = [
    "sexo",
    "obesidad"
]


# ==========================================================================
# AUXILIARY FUNCTIONS
# ==========================================================================

def convert_to_numeric(series: pd.Series) -> pd.Series:
    """
    Convert a pandas Series to numeric format.

    Supports decimal commas by replacing ',' with '.'.
    """

    if pd.api.types.is_numeric_dtype(series):
        return series

    s = series.astype(str).str.strip()
    s = s.str.replace(",", ".", regex=False)

    return pd.to_numeric(s, errors="coerce")


# ==========================================================================
# STATISTICAL TESTS
# ==========================================================================

def welch_t_test(x0: np.ndarray, x1: np.ndarray) -> float:
    """
    Perform Welch's t-test.
    """

    x0 = x0[~np.isnan(x0)]
    x1 = x1[~np.isnan(x1)]

    if x0.size < 2 or x1.size < 2:
        return np.nan

    return float(
        stats.ttest_ind(
            x0,
            x1,
            equal_var=False,
            nan_policy="omit"
        ).pvalue
    )


def mann_whitney_u_test(x0: np.ndarray, x1: np.ndarray) -> float:
    """
    Perform Mann–Whitney U test.
    """

    x0 = x0[~np.isnan(x0)]
    x1 = x1[~np.isnan(x1)]

    if x0.size < 1 or x1.size < 1:
        return np.nan

    return float(
        stats.mannwhitneyu(
            x0,
            x1,
            alternative="two-sided",
            method="asymptotic"
        ).pvalue
    )


def chi_square_or_fisher(
    contingency_table: np.ndarray
) -> tuple[str, float, np.ndarray]:
    """
    Perform Chi-square or Fisher's exact test.

    Fisher's exact test is only applied to 2x2 tables
    when expected frequencies are < 5.
    """

    chi2, p_chi, dof, expected = stats.chi2_contingency(
        contingency_table,
        correction=False
    )

    if (
        contingency_table.shape == (2, 2)
        and (expected < 5).any()
    ):

        _, p_fisher = stats.fisher_exact(
            contingency_table,
            alternative="two-sided"
        )

        return "Fisher_Exact_Test", float(p_fisher), expected

    else:

        return "Chi_Square_Test", float(p_chi), expected


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
        help="Input CSV/TSV file containing T2D status."
    )

    parser.add_argument(
        "--sep",
        default="\t",
        help="File separator. Default = tab-separated."
    )

    parser.add_argument(
        "--outdir",
        default="statistical_tests_output",
        help="Output directory."
    )

    parser.add_argument(
        "--group-col",
        default="DM2",
        help="Grouping column (default = DM2)."
    )

    args = parser.parse_args()

    # ----------------------------------------------------------------------
    # CREATE OUTPUT DIRECTORY
    # ----------------------------------------------------------------------

    os.makedirs(args.outdir, exist_ok=True)

    # ----------------------------------------------------------------------
    # LOAD INPUT DATA
    # ----------------------------------------------------------------------

    df = pd.read_csv(
        args.infile,
        sep=args.sep,
        dtype=str
    )

    # ----------------------------------------------------------------------
    # CHECK GROUP COLUMN
    # ----------------------------------------------------------------------

    if args.group_col not in df.columns:

        raise ValueError(
            f"Grouping column '{args.group_col}' not found."
        )

    df[args.group_col] = pd.to_numeric(
        df[args.group_col],
        errors="coerce"
    )

    # ----------------------------------------------------------------------
    # SPLIT BY T2D STATUS
    # ----------------------------------------------------------------------

    controls = df[df[args.group_col] == 0].copy()
    cases = df[df[args.group_col] == 1].copy()

    # ----------------------------------------------------------------------
    # RESULTS STORAGE
    # ----------------------------------------------------------------------

    statistical_results = []
    contingency_results = []

    # ======================================================================
    # MEAN ± SD VARIABLES → WELCH'S T-TEST
    # ======================================================================

    for variable in MEAN_SD_VARIABLES:

        if variable not in df.columns:
            continue

        x0 = convert_to_numeric(
            controls[variable]
        ).to_numpy(dtype=float)

        x1 = convert_to_numeric(
            cases[variable]
        ).to_numpy(dtype=float)

        n_controls = int(np.sum(~np.isnan(x0)))
        n_cases = int(np.sum(~np.isnan(x1)))

        p_value = welch_t_test(x0, x1)

        statistical_results.append({

            "Variable": variable,
            "Variable_Type": "Quantitative_Mean_SD",

            "Statistical_Test": "Welch_t_test",

            "P_Value": p_value,

            "N_Controls": n_controls,
            "N_Cases": n_cases
        })

    # ======================================================================
    # MEDIAN [IQR] VARIABLES → MANN-WHITNEY U
    # ======================================================================

    for variable in MEDIAN_IQR_VARIABLES:

        if variable not in df.columns:
            continue

        x0 = convert_to_numeric(
            controls[variable]
        ).to_numpy(dtype=float)

        x1 = convert_to_numeric(
            cases[variable]
        ).to_numpy(dtype=float)

        n_controls = int(np.sum(~np.isnan(x0)))
        n_cases = int(np.sum(~np.isnan(x1)))

        p_value = mann_whitney_u_test(x0, x1)

        statistical_results.append({

            "Variable": variable,
            "Variable_Type": "Quantitative_Median_IQR",

            "Statistical_Test": "Mann_Whitney_U_Test",

            "P_Value": p_value,

            "N_Controls": n_controls,
            "N_Cases": n_cases
        })

    # ======================================================================
    # CATEGORICAL VARIABLES
    # ======================================================================

    for variable in CATEGORICAL_VARIABLES:

        if variable not in df.columns:
            continue

        y0 = pd.to_numeric(
            controls[variable],
            errors="coerce"
        )

        y1 = pd.to_numeric(
            cases[variable],
            errors="coerce"
        )

        categories = sorted(
            pd.concat([y0, y1], axis=0)
            .dropna()
            .unique()
            .tolist()
        )

        if len(categories) == 0:
            continue

        # ------------------------------------------------------------------
        # BUILD CONTINGENCY TABLE
        # ------------------------------------------------------------------

        contingency_table = []

        for category in categories:

            contingency_table.append([
                int((y0 == category).sum()),
                int((y1 == category).sum())
            ])

        contingency_table = np.array(
            contingency_table,
            dtype=int
        )

        # ------------------------------------------------------------------
        # SELECT APPROPRIATE TEST
        # ------------------------------------------------------------------

        test_name, p_value, expected = chi_square_or_fisher(
            contingency_table
        )

        statistical_results.append({

            "Variable": variable,

            "Variable_Type": (
                "Binary_Categorical"
                if len(categories) == 2
                else "Multilevel_Categorical"
            ),

            "Statistical_Test": test_name,

            "P_Value": p_value,

            "N_Controls": int(y0.notna().sum()),
            "N_Cases": int(y1.notna().sum()),

            "Categories": ",".join(map(str, categories))
        })

        # ------------------------------------------------------------------
        # STORE CONTINGENCY TABLE
        # ------------------------------------------------------------------

        contingency_results.append({

            "Variable": variable,

            "Categories": ",".join(
                map(str, categories)
            ),

            "Contingency_Table": ";".join([
                f"{categories[i]}:"
                f"{contingency_table[i,0]},"
                f"{contingency_table[i,1]}"
                for i in range(len(categories))
            ])
        })

    # ----------------------------------------------------------------------
    # SAVE STATISTICAL RESULTS
    # ----------------------------------------------------------------------

    statistical_results_df = (
        pd.DataFrame(statistical_results)
        .sort_values(["Variable_Type", "Variable"])
    )

    statistical_output = os.path.join(
        args.outdir,
        "statistical_tests_results.csv"
    )

    statistical_results_df.to_csv(
        statistical_output,
        index=False
    )

    # ----------------------------------------------------------------------
    # SAVE CONTINGENCY TABLES
    # ----------------------------------------------------------------------

    contingency_results_df = pd.DataFrame(
        contingency_results
    )

    contingency_output = os.path.join(
        args.outdir,
        "contingency_tables.csv"
    )

    contingency_results_df.to_csv(
        contingency_output,
        index=False
    )

    # ----------------------------------------------------------------------
    # FINAL REPORT
    # ----------------------------------------------------------------------

    print("\n====================================================")
    print("Statistical tests successfully completed")
    print("====================================================")

    print(f"\nStatistical results:")
    print(f"  {statistical_output}")

    print(f"\nContingency tables:")
    print(f"  {contingency_output}")

    print("\nSTATISTICAL TESTS APPLIED:")
    print("- Mean ± SD variables → Welch's t-test")
    print("- Median [IQR] variables → Mann–Whitney U test")
    print("- Categorical variables → Chi-square test")
    print("- Fisher's exact test applied automatically")
    print("  when expected frequencies < 5 in 2x2 tables")


# ==========================================================================
# SCRIPT EXECUTION
# ==========================================================================

if __name__ == "__main__":
    main()
