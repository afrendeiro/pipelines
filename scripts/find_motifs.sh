#! /bin/bash

REFDIR="~/data/human/E2F7/"
OUTDIR="~/output/human/E2F7/"
GENOMEREF="/home/s3/afr/data/human/annotation/hg19/hg19.fa"
GENOMESIZES="/home/s3/afr/data/human/annotation/hg19/hg19_ChrSizes.tsv"
MOTIFDB="/home/s3/afr/bin/meme_4.9.1/db/motif_databases/JASPAR_CORE_2009_vertebrates.meme"

# get 500bp region around peak center
bedtools slop -i ${OUTDIR}/E2F7_Hela.peaks.bed -g $GENOMESIZES -b 200 |

# get fasta sequences from bed files
bedtools getfasta -fi $GENOMEREF -bed - -fo ${REFDIR}/E2F7_Hela.peaks.fa

# small curated database
meme-chip -meme-p 12 -db $MOTIFDB  -oc ${OUTDIR}/meme-chip ${REFDIR}/E2F7_Hela.peaks.fa

# for urochordates
JASPAR_CORE_2009_urochordates.meme