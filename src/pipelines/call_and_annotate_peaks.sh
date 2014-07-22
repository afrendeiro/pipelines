#!/bin/bash

## CALL PEAKS

MAPPED=/sysdev/s3/share/data/oikopleura/chip-seq/mapped
PEAKS=/sysdev/s3/share/data/oikopleura/chip-seq/peaks

function call_peaks {
    STAGE=$1
    MAPPED=$2
    PEAKS=$3

    for SAMPLE in ${STAGE}_{'E2F1','E2F7','POL2'}_R1.bwa
    do
        if [ ! -s ${MAPPED}/${SAMPLE}.peaks.tsv ]; then
            echo `date +%Y.%m.%d-%H:%M:%S` - "Finding peaks for $SAMPLE..." >> log.txt
            INPUT=${STAGE}_IN_R1.bwa
            # Run peakzilla
            peakzilla.py ${MAPPED}/${SAMPLE}.ext.bed ${MAPPED}/${INPUT}.bed > ${PEAKS}/${SAMPLE}.peaks.tsv
            # Make bed file
            cut -f 1,2,3 ${PEAKS}/${SAMPLE}.peaks.tsv | tail -n +2 > ${PEAKS}/${SAMPLE}.peaks.bed
        fi
    done
}

export -f call_peaks

parallel call_peaks ::: {'TB','D2','D6-F','D6-M'} ::: $MAPPED ::: $PEAKS

# ANNOTATE PEAKS

PEAKS=/sysdev/s3/share/data/oikopleura/chip-seq/peaks

CHRSIZES=~/data/oikopleura/assembly/Oikopleura_reference_chrSizes.tsv
GENOMELOCS=~/data/oikopleura/annotation/Oikopleura_gene_models.genomic_locations.bed
TSSs=~/data/oikopleura/annotation/Oikopleura_gene_models_reference.TSSs.bed

function annotate {
    SAMPLE=$1
    SAMPLE_SHUFFLED=${SAMPLE}.shuffled
    PEAKS=$2
    CHRSIZES=$3
    GENOMELOCS=$4
    TSSs=$5

    # generate shuffled peaks
    if [ -s ${PEAKS}/${SAMPLE}.bed ] && [ ! -s ${PEAKS}/${SAMPLE_SHUFFLED}.bed ]; then
        echo `date +%Y.%m.%d-%H:%M:%S` - "Generating random peaks for $SAMPLE..." >> log.txt
        bedtools shuffle -i ${PEAKS}/${SAMPLE}.tsv -g $CHRSIZES -seed 2014 | bedtools sort > ${PEAKS}/${SAMPLE_SHUFFLED}.tsv
    fi

    # Annotate with genomic features and nearest TSS
    ## reports an entry for all features (even not intercepting, 0)
    if [ -s ${PEAKS}/${SAMPLE}.tsv ] && [ -s ${PEAKS}/${SAMPLE_SHUFFLED}.tsv ] && [ ! -s ${PEAKS}/${SAMPLE}.annotation.tsv ] && [ ! -s ${PEAKS}/${SAMPLE_SHUFFLED}.annotation.tsv ]; then
        echo `date +%Y.%m.%d-%H:%M:%S` - "Annotating $SAMPLE with annotation..." >> log.txt
        bedtools intersect -wao -a ${PEAKS}/${SAMPLE}.tsv -b $GENOMELOCS | bedtools closest -d -t all -a - -b $TSSs > ${PEAKS}/${SAMPLE}.annotation.tsv

        echo `date +%Y.%m.%d-%H:%M:%S` - "Annotating $SAMPLE_SHUFFLED with annotation..." >> log.txt
        bedtools intersect -wao -a ${PEAKS}/${SAMPLE_SHUFFLED}.tsv -b $GENOMELOCS | bedtools closest -d -t all -a - -b $TSSs > ${PEAKS}/${SAMPLE_SHUFFLED}.annotation.tsv
    fi

    # Annotate D6 samples with chromatin states
    if [[ $SAMPLE == *D6* ]]; then
        # choose appropriate chromatin file 
        if [[ $SAMPLE == *D6-F* ]]
            then
            CHROMATIN=~/data/oikopleura/annotation/chromHMM/Female_15_segments.bed
            else
            CHROMATIN=~/data/oikopleura/annotation/chromHMM/Male_15_segments.bed
        fi
        R=$RANDOM

        echo `date +%Y.%m.%d-%H:%M:%S` - "Annotating $SAMPLE with chromatin states..." >> log.txt
        bedtools intersect -wao -f 0.3 -a ${PEAKS}/${SAMPLE}.annotation.tsv -b $CHROMATIN > TMP_${R}
        mv TMP_${R} ${PEAKS}/${SAMPLE}.annotation.tsv

        echo `date +%Y.%m.%d-%H:%M:%S` - "Annotating $SAMPLE_SHUFFLED with chromatin states..." >> log.txt
        bedtools intersect -wao -f 0.3 -a ${PEAKS}/${SAMPLE_SHUFFLED}.annotation.tsv -b $CHROMATIN > TMP_${R}
        mv TMP_${R} ${PEAKS}/${SAMPLE_SHUFFLED}.annotation.tsv
    fi
}

export -f annotate

parallel annotate ::: {'TB','D2','D6-F','D6-M'}_{'E2F1','E2F7','POL2'}_R1.bwa.peaks ::: $PEAKS ::: $CHRSIZES ::: $GENOMELOCS ::: $TSSs
