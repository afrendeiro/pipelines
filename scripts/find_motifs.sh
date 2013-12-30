#! /bin/bash

REFDIR="~/data/human/E2F7/"
OUTDIR="~/output/human/E2F7/"

# get fasta sequences from bed files
bedtools getfasta -fi ~/data/human/annotation/hg19/hg19.fa -bed ${REFDIR}/E2F7_Hela.peaks.bed -fo ${REFDIR}/E2F7_Hela.peaks.fa

# small curated database
meme-chip -meme-p 12 -db ~/bin/meme_4.9.1/db/motif_databases/JASPAR_CORE_2009_vertebrates.meme  -oc ${OUTDIR}/meme-chip ${REFDIR}/E2F7_Hela.peaks.fa

# for urochordates
JASPAR_CORE_2009_urochordates.meme