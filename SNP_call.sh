#!/usr/bin/env bash
set -euo pipefail

REF="CHO_EColi_Dul-6plasm.fa"
REGIONS="p1.1-Tr2-Dul,p1.2-GS-Dul"

BAMS=(
  R12_Ref1.bam
  MCB_Ref1.bam
  WCB_Ref1.bam
  EoP_Ref1.bam
)

MIN_ALT=1
MIN_DP=5

OUTDIR="vcf3"
mkdir -p "$OUTDIR"

OUT_CSV="${OUTDIR}/all_samples_minimal.csv"
echo "Sample,CHROM,POS,REF,ALT,QUAL,FORMAT" > "$OUT_CSV"

### Variant calling and export

for BAM in "${BAMS[@]}"; do
    SAMPLE="${BAM%.bam}"
    VCF="${OUTDIR}/${SAMPLE}.vcf"

    samtools index "$BAM" || true

    bcftools mpileup \
        -f "$REF" \
        -r "$REGIONS" \
        -a FORMAT/AD,FORMAT/DP \
        -Q 10 -q 10 \
        -d 5000 \
        "$BAM" \
    | bcftools call -m --ploidy 1 \
    | bcftools filter -i "FORMAT/DP>=${MIN_DP} && FORMAT/AD[0:1]>=${MIN_ALT}" \
    > "$VCF"

    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT:%PL:%DP:%AD]\n' \
        "$VCF" \
    | awk -v sample="$SAMPLE" -F'\t' \
        '{print sample "," $1 "," $2 "," $3 "," $4 "," $5 "," $6}' \
    >> "$OUT_CSV"
done
