#!/bin/bash

MAPPED=/sysdev/s3/share/data/oikopleura/chip-seq/mapped
NORMALIZED=/sysdev/s3/share/data/oikopleura/chip-seq/normalized

WINDOWS=~/data/oikopleura/assembly/50bp_25bp_sliding_windows.bed
CHRSIZES=~/data/oikopleura/assembly/Oikopleura_reference_chrSizes.tsv

# Make windows genome-wide
bedtools makewindows -g $CHRSIZES -w 50 -s 25 > $WINDOWS

# Calculate coverage in those windows for each sample
for SAMPLE in {'TB','D2','D6-F','D6-M'}_{'IN','E2F1','E2F7','POL2'}_R1.bwa
do
    if [ -s ${MAPPED}/${SAMPLE}.coverage.bed ]; then
        echo `date +%Y.%m.%d-%H:%M:%S` - "Coverage for $SAMPLE already exist: ${SAMPLE}.bwa.bam, skipping.." >> log.txt
    else
        echo `date +%Y.%m.%d-%H:%M:%S` - "Calculating coverage for $SAMPLE..." >> log.txt
        bedtools coverage -a ${MAPPED}/${SAMPLE}.ext.bed -b $WINDOWS > ${MAPPED}/${SAMPLE}.coverage.bed
    fi
done

# Concatenate all files, normalize and export bedgraphs
R CMD concatenate_and_normalize.R

# Make BigWig tracks
for SAMPLE in {'TB','D2','D6-F','D6-M'}_{'E2F1','E2F7','POL2'}_R1.bwa
do
    if [ -s ${NORMALIZED}/${SAMPLE}.norm.bw ]; then
        echo `date +%Y.%m.%d-%H:%M:%S` - "BigWig for $SAMPLE already exist: ${SAMPLE}.norm.bw, skipping.." >> log.txt
    else
        echo `date +%Y.%m.%d-%H:%M:%S` - "Making BigWig for normalized $SAMPLE.norm..." >> log.txt
        TMP=tmp_$RANDOM
        bedtools sort -i ${NORMALIZED}/${SAMPLE}.norm.bedgraph > $TMP
        mv $TMP ${NORMALIZED}/${SAMPLE}.norm.bedgraph
        bedGraphToBigWig ${NORMALIZED}/${SAMPLE}.norm.bedgraph $CHRSIZES ${NORMALIZED}/${SAMPLE}.norm.bw
    fi
done