#!/bin/bash
#SBATCH -p comet
#SBATCH --job-name LabMCB_20250306
#SBATCH -c 128

REF="/mnt/storage/lab13/Ref_2/Ref2.fa"
READ="/mnt/storage/lab13/R3.fastq.gz"
OUT2="/mnt/storage/lab13/Ref_2/R3_Ref2.bam"

/mnt/storage/lab13/minimap2-2.28_x64-linux/minimap2 -a -t 128 $REF $READ | /mnt/storage/lab13/samtools-1.20/samtools view -Sbh -F4 | /mnt/storage/lab13/samtools-1.20/samtools sort -@ 128 > $OUT2

/mnt/storage/lab13/samtools-1.20/samtools index $OUT2

/mnt/storage/lab13/samtools-1.20/samtools flagstat $OUT2 > "$OUT2"_stat.txt

READS1="/mnt/storage/lab13/Dul_PE/250123_HSGA.Vorobev_DNA_2024.DUL73.read1.fastq.gz"
READS2="/mnt/storage/lab13/Dul_PE/250123_HSGA.Vorobev_DNA_2024.DUL73.read2.fastq.gz"
OUT1="/mnt/storage/lab13/Ref_2/R12_Ref2.bam"

/mnt/storage/lab13/bwa/bwa index $REF

/mnt/storage/lab13/bwa/bwa mem -t 128 -R "@RG\tID:id1\tSM:sample1\tLB:lib1" $REF $READS1 $READS2 | /mnt/storage/lab13/samtools-1.20/samtools view -Sbh -F4 | /mnt/storage/lab13/samtools-1.20/samtools sort -@ 128 > $OUT1
/mnt/storage/lab13/samtools-1.20/samtools index $OUT1
/mnt/storage/lab13/samtools-1.20/samtools flagstat $OUT1 > "$OUT1"_stat.txt

READS3="/mnt/storage/lab13/4BGD/ZDEISCR2-GED-G007WGS-MG-350129749-2-69-F.fq.gz"
READS4="/mnt/storage/lab13/4BGD/ZDEISCR2-GED-G007WGS-MG-350129749-2-69-R.fq.gz"
OUT3="/mnt/storage/lab13/Ref_2/4BGD_Ref2.bam"

/mnt/storage/lab13/bwa/bwa mem -t 128 -R "@RG\tID:id1\tSM:sample1\tLB:lib1" $REF $READS3 $READS4 | /mnt/storage/lab13/samtools-1.20/samtools view -Sbh -F4 | /mnt/storage/lab13/samtools-1.20/samtools sort -@ 128 > $OUT3
/mnt/storage/lab13/samtools-1.20/samtools index $OUT3
/mnt/storage/lab13/samtools-1.20/samtools flagstat $OUT3 > "$OUT3"_stat.txt
