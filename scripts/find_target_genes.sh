#! /bin/bash

# Finds target genes based on TFBS based on proximity to nearest TSS

REFDIR="/home/s3/afr/data/human/E2F7"
OUTDIR="/home/s3/afr/output/human/E2F7"

# Extract peak summit position 
# Find closest gene
# measure distance to TSS
# filter where distance is < X
tail -n +2 ${OUTDIR}/E2F7_Hela.peaks.tsv | awk '{print $1,$5,($5+1)}' | sed '{:q;N;s/ /\t/g;t q}' | \
bedtools closest -d -a - -b ~/data/human/annotation/hg19/genomic_location/TSS.bed | \
cut -f 1,2,3,7,10 > ${OUTDIR}/E2F7_Hela.targets.bed

cut -f 4 ${OUTDIR}/E2F7_Hela.targets.bed | sed 's/\([^_]\_[^_]*\).*/\1/' > ${OUTDIR}/E2F7_Hela.targets.txt



# Find out if gene is part of Operon
