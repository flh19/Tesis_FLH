#!/bin/bash

# ============================================================
# SCRIPT: extracción de métricas QC del VCF
# ============================================================
#
# Métricas:
#   - DP  (Depth per genotype)
#   - GQ  (Genotype Quality)
#   - MQ  (Mapping Quality)
#
# Se generan:
#   1. Valores individuales por genotipo
#   2. Valores globales en una columna
#   3. Media por variante
#
# Uso:
#   bash qc_metrics.sh BEFORE
#   bash qc_metrics.sh AFTER
#
# Ejemplo:
#   bash qc_metrics.sh before
#   bash qc_metrics.sh after
#
# ============================================================

set -euo pipefail

# ============================================================
# INPUTS
# ============================================================

STATUS=$1

VCF="VCexomes_*.vcf.bgz"

OUTDIR="QC_${STATUS}"

mkdir -p ${OUTDIR}

echo "======================================="
echo "Running QC extraction: ${STATUS}"
echo "======================================="

# ============================================================
# 1. DP (Depth)
# ============================================================

echo "Extracting DP..."

bcftools query \
    -f '[%DP\t]\n' \
    ${VCF} \
    > ${OUTDIR}/DP_genotypes_${STATUS}.txt

# ------------------------------------------------------------
# Convertir a una sola columna
# ------------------------------------------------------------

tr '\t' '\n' \
    < ${OUTDIR}/DP_genotypes_${STATUS}.txt \
    | grep -v '\.' \
    > ${OUTDIR}/DP_all_values_${STATUS}.txt

# ------------------------------------------------------------
# Media DP por variante
# ------------------------------------------------------------

awk '
{
    sum=0
    n=0

    for(i=1;i<=NF;i++){

        if($i!="."){

            sum += $i
            n++
        }
    }

    if(n>0){

        print sum/n

    } else {

        print "."
    }
}
' ${OUTDIR}/DP_genotypes_${STATUS}.txt \
> ${OUTDIR}/DP_mean_per_variant_${STATUS}.txt

# ============================================================
# 2. GQ (Genotype Quality)
# ============================================================

echo "Extracting GQ..."

bcftools query \
    -f '[%GQ\t]\n' \
    ${VCF} \
    > ${OUTDIR}/GQ_genotypes_${STATUS}.txt

# ------------------------------------------------------------
# Convertir a una sola columna
# ------------------------------------------------------------

tr '\t' '\n' \
    < ${OUTDIR}/GQ_genotypes_${STATUS}.txt \
    | grep -v '\.' \
    > ${OUTDIR}/GQ_all_values_${STATUS}.txt

# ------------------------------------------------------------
# Media GQ por variante
# ------------------------------------------------------------

awk '
{
    sum=0
    n=0

    for(i=1;i<=NF;i++){

        if($i!="."){

            sum += $i
            n++
        }
    }

    if(n>0){

        print sum/n

    } else {

        print "."
    }
}
' ${OUTDIR}/GQ_genotypes_${STATUS}.txt \
> ${OUTDIR}/GQ_mean_per_variant_${STATUS}.txt

# ============================================================
# 3. MQ (Mapping Quality)
# ============================================================

echo "Extracting MQ..."

bcftools query \
    -f '%MQ\n' \
    ${VCF} \
    > ${OUTDIR}/MQ_variant_${STATUS}.txt

grep -v '\.' \
    ${OUTDIR}/MQ_variant_${STATUS}.txt \
    > ${OUTDIR}/MQ_variant_clean_${STATUS}.txt

# ============================================================
# RESUMEN
# ============================================================

echo ""
echo "======================================="
echo "QC extraction completed"
echo "======================================="

echo ""
echo "Generated files:"
echo ""

ls -lh ${OUTDIR}

echo ""
echo "======================================="
