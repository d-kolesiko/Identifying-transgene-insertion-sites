#!/usr/bin/env bash
set -euo pipefail

BAM="R3_Ref2.bam"
PLASMID="p1.2-GS-Dul"

MIN_MAPQ=10
MIN_LEN=3000

OUTDIR="p12_res"
mkdir -p "$OUTDIR"

READ_LIST="${OUTDIR}/${PLASMID}_reads.lst"
FINAL_CSV="${OUTDIR}/${PLASMID}_intervals.csv"

### Stage 1: reads aligned to plasmid

samtools view "$BAM" "$PLASMID" \
| awk '{print $1}' \
| sort -u \
> "$READ_LIST"

echo "Stage 1: $(wc -l < "$READ_LIST") reads with plasmid alignment"

### Stage 2: filter alignments and collect intervals

samtools view "$BAM" \
| awk -v plasmid="$PLASMID" \
      -v min_mapq="$MIN_MAPQ" \
      -v min_len="$MIN_LEN" \
      -v id_file="$READ_LIST" '
BEGIN {
    FS = "\t"
    while ((getline < id_file) > 0) keep[$1] = 1
}
{
    read_id = $1
    if (!(read_id in keep)) next

    ref   = $3
    start = $4
    mapq  = $5
    seq   = $10

    if (mapq <= min_mapq) next
    len = length(seq)
    if (len < min_len) next

    end = start + len - 1
    key = read_id SUBSEP ref SUBSEP start SUBSEP end
    if (seen[key]++) next

    if (ref == plasmid) {
        pi = ++plasmid_n[read_id]
        plasmid_ref[read_id,pi]   = ref
        plasmid_start[read_id,pi] = start
        plasmid_end[read_id,pi]   = end
    } else {
        gi = ++genome_n[read_id]
        genome_ref[read_id,gi]   = ref
        genome_start[read_id,gi] = start
        genome_end[read_id,gi]   = end
    }

    has_read[read_id] = 1
}
END {
    valid = 0
    for (r in has_read)
        if ((r in plasmid_n) && (r in genome_n))
            valid++

    print "Stage 2: " valid " reads after MAPQ and length filters" > "/dev/stderr"

    out = "'"$FINAL_CSV"'"
    print "read_id,interval,length" > out

    for (r in has_read) {
        if (!(r in plasmid_n) || !(r in genome_n)) continue

        for (i = 1; i <= plasmid_n[r]; i++) {
            s = plasmid_start[r,i]
            e = plasmid_end[r,i]
            print r "," plasmid_ref[r,i] ":" s "-" e "," e - s + 1 >> out
        }

        for (i = 1; i <= genome_n[r]; i++) {
            s = genome_start[r,i]
            e = genome_end[r,i]
            print r "," genome_ref[r,i] ":" s "-" e "," e - s + 1 >> out
        }
    }
}
'
