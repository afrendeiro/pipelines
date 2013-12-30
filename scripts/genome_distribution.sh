#! /bin/bash

# Genomic Regions:
	#TSS (-1kb to +100bp)
	#TTS (-100 bp to +1kb)
	#CDS Exons
	#5' UTR
	#3' UTR
	#**CpG Islands
	#**Repeats
	#Introns (1st intron?)
	#Intergenic

REFDIR="~/data/human/E2F7/"
OUTDIR="~/output/human/E2F7/"

# Intersect 
# reports an entry for all features (even not intercepting, 0)
bedtools intersect -wao -a ${REFDIR}/E2F7_Hela.peaks.bed -b ~/data/human/annotation/hg19/hg19.location.bed | cut -f 1,2,3,7,10 > ${OUTDIR}/E2F7_Hela.peaks.location.bed

# measure distance to TSS
bedtools closest -d -a ${REFDIR}/E2F7_Hela.peaks.bed -b ~/data/human/annotation/hg19/genomic_location/TSSs.bed | cut -f 1,2,3,7,10 > ${OUTDIR}/E2F7_Hela.peaks.TSS.bed