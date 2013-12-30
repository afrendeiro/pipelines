#!/bin/bash

REFDIR="~/data/human/E2F7"
OUTDIR="~/output/human/E2F7"

peakzilla.py ${REFDIR}/GSM810993_E2F7_Hela.081213.bwa.bed GSM810994_Input_Hela.081213.bwa.bed > ${OUTDIR}/E2F7_Hela.peaks.tsv
cut -f 1,2,3 E2F7_Hela.peaks.tsv | tail -n +2 > ${OUTDIR}/E2F7_Hela.peaks.bed

# compare peak sets


