#!/bin/bash
# Run meme on samples

MAPPED=/sysdev/s3/share/data/oikopleura/chip-seq/mapped
PEAKS=/sysdev/s3/share/data/oikopleura/chip-seq/peaks
MOTIFS=/sysdev/s3/share/data/oikopleura/chip-seq/motifs

GENOMEREF=~/data/oikopleura/assembly/Oikopleura_reference_masked_v3.0.fa
CHRSIZES=~/data/oikopleura/assembly/Oikopleura_reference_chrSizes.tsv
MEMEDB=~/memeDBs/JASPAR_CORE_2014_vertebrates.meme

for SAMPLE in {'TB','D2','D6-F','D6-M'}_{'E2F1','E2F7'}_R1.bwa
do
	# get peaks sumits in bed format
	echo `date +%Y.%m.%d-%H:%M:%S` - "Getting summits of $SAMPLE peaks..." >> log.txt
	tail -n +2 ${PEAKS}/$SAMPLE.peaks.tsv | awk -F OFS="\t" '{print $1,$5,($5+1)}' > ${PEAKS}/$SAMPLE.peaks.summits.bed

	# get 500bp region around peak summit and get fasta sequences
	echo `date +%Y.%m.%d-%H:%M:%S` - "Getting 500bp around summit of $SAMPLE peaks..." >> log.txt
	bedtools slop -i ${PEAKS}/$SAMPLE.peaks.summits.bed -g $CHRSIZES -b 250 | bedtools getfasta -fi $GENOMEREF -bed - -fo ${PEAKS}/$SAMPLE.peaks.summits.fa

	# make output directory
	echo `date +%Y.%m.%d-%H:%M:%S` - "Making directory ${MOTIFS}/${SAMPLE}-meme..." >> log.txt
	mkdir ${MOTIFS}/${SAMPLE}-meme

	# small curated database
	echo `date +%Y.%m.%d-%H:%M:%S` - "Running meme-chip on $SAMPLE..." >> log.txt
	meme-chip -meme-p 24 -db $MEMEDB -oc ${MOTIFS}/${SAMPLE}-meme ${PEAKS}/$SAMPLE.peaks.summits.fa
done
