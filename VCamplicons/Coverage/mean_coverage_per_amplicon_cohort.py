#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Per-amplicon coverage summary for a target gene.

Description:
This script processes sequencing coverage files from multiple cohorts
and calculates per-amplicon coverage statistics for a user-defined gene.

The workflow is fully generalizable to any target gene by modifying
the GEN_OBJETIVO variable.
"""

# =========================================================
# Imports
# =========================================================

import os

import pandas as pd


# =========================================================
# Global configuration
# =========================================================

RUTAS = {

    "Di@bet.es": "../Diabet_study/cov_files/",

    "Hortega": "../Hortega_study/cov_files/"
}

# ---------------------------------------------------------
# Target gene
# ---------------------------------------------------------
# Modify this variable to analyze another gene

GEN_OBJETIVO = "NAT1"

# ---------------------------------------------------------
# Output file
# ---------------------------------------------------------

OUTFILE = (
    f"nivel3_cobertura_media_por_amplicon_"
    f"{GEN_OBJETIVO}.tsv"
)


# =========================================================
# Coverage processing functions
# =========================================================

def cobertura_media_por_amplicon(
    ruta,
    poblacion
):
    """
    Process all .cov files from a cohort directory and
    calculate per-amplicon coverage statistics for the
    target gene.

    Parameters
    ----------
    ruta : str
        Path to the directory containing .cov files.

    poblacion : str
        Cohort name.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing:
        - amplicon identifier
        - coverage summary statistics
        - cohort name
    """

    resultados = []

    # -----------------------------------------------------
    # Iterate through coverage files
    # -----------------------------------------------------

    for fichero in os.listdir(ruta):

        if not fichero.endswith(".cov"):
            continue

        # -------------------------------------------------
        # Read coverage file
        # -------------------------------------------------

        df = pd.read_csv(
            os.path.join(ruta, fichero),

            sep=r"\s+",

            header=None,

            names=[
                "muestra",
                "amplicon",
                "posicion",
                "cobertura"
            ]
        )

        # -------------------------------------------------
        # Ensure numeric coverage values
        # -------------------------------------------------

        df["cobertura"] = pd.to_numeric(
            df["cobertura"],
            errors="coerce"
        ).fillna(0)

        # -------------------------------------------------
        # Extract gene name from amplicon identifier
        # -------------------------------------------------

        df["gen"] = (
            df["amplicon"]
            .str.split("-")
            .str[0]
        )

        # -------------------------------------------------
        # Filter target gene
        # -------------------------------------------------

        df = df[
            df["gen"] == GEN_OBJETIVO
        ]

        # -------------------------------------------------
        # Per-amplicon coverage statistics
        # -------------------------------------------------

        resumen = (

            df.groupby("amplicon")["cobertura"]

            .agg(

                media="mean",

                sd="std",

                mediana="median",

                p25=lambda x:
                x.quantile(0.25),

                p27=lambda x:
                x.quantile(0.27)
            )

            .reset_index()
        )

        # -------------------------------------------------
        # Add cohort metadata
        # -------------------------------------------------

        resumen["poblacion"] = poblacion

        resultados.append(resumen)

    return pd.concat(
        resultados,
        ignore_index=True
    )


# =========================================================
# Coverage calculation
# =========================================================

lista = []

for poblacion, ruta in RUTAS.items():

    print(
        f"Processing {poblacion}..."
    )

    lista.append(
        cobertura_media_por_amplicon(
            ruta,
            poblacion
        )
    )

# ---------------------------------------------------------
# Merge cohort results
# ---------------------------------------------------------

df_total = pd.concat(
    lista,
    ignore_index=True
)

# ---------------------------------------------------------
# Global summary statistics
# ---------------------------------------------------------

resumen_final = (

    df_total

    .groupby(
        ["poblacion", "amplicon"]
    )

    .agg(

        cobertura_media=(
            "media",
            "mean"
        ),

        cobertura_sd=(
            "sd",
            "mean"
        ),

        cobertura_mediana=(
            "mediana",
            "mean"
        ),

        cobertura_p25=(
            "p25",
            "mean"
        ),

        cobertura_p27=(
            "p27",
            "mean"
        )
    )

    .reset_index()
)


# =========================================================
# Export results
# =========================================================

resumen_final.to_csv(
    OUTFILE,

    sep="\t",

    index=False
)


# =========================================================
# Final message
# =========================================================

print(
    f"Output file generated: "
    f"{OUTFILE}"
)
