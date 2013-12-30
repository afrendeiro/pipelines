#!/bin/bash

GENOMEREF="data/human/annotation/hg19/hg19.fa"
GENOMESIZES="data/human/annotation/hg19/hg19_ChrSizes.tsv"
ROOTDIR="/home/s3/afr/data/human/E2F7"
OUT_DIR="/home/s3/afr/data/human/E2F7/q30"

run_bwa_pipeline () {
    local DIR=$1
    local INFILE=$2
    local EXPERIMENT_NAME=$3
    local OUTDIR=$4

    if [ -s ${OUTDIR}/${EXPERIMENT_NAME}.bwa.bam ]; then
        echo "Results for $EXPERIMENT_NAME already exist: ${OUTDIR}/${EXPERIMENT_NAME}.bwa.bam, skipping.."
    else
        echo "Running bwa alignment for $EXPERIMENT_NAME.."  #output to stdout
    	echo "FILE:$DIR/$INFILE NAME:$EXPERIMENT_NAME OUTDIR:$OUTDIR .."

        #align to genome, output .bam file
        # filters for quality >=q30
    	bwa aln -t 8 $GENOMEREF ${DIR}/${INFILE}.fastq | bwa samse $GENOMEREF - ${DIR}/${INFILE} | samtools view -b -q 30 -u -T $GENOMEREF - | samtools sort - ${OUTDIR}/${EXPERIMENT_NAME}.bwa.bam
    	
        # create index
        samtools index ${OUTDIR}/${EXPERIMENT_NAME}.bwa.bam
        mv ${OUTDIR}/${EXPERIMENT_NAME}.bwa.bam.bai ${OUTDIR}/${EXPERIMENT_NAME}.bwa.bai

        # count reads, export to file
        samtools view -c ${OUTDIR}/${EXPERIMENT_NAME}.bwa.bam | awk -v var=$EXPERIMENT_NAME '$2 = $2 FS var' > ${OUTDIR}/${INFILE}.readcount.txt

        # convert to bed
        bamToBed -i ${OUTDIR}/${EXPERIMENT_NAME}.bwa.bam ${OUTDIR}/${EXPERIMENT_NAME}.bwa.bed

        # extend reads
        slopBed -g $GENOMESIZES -l 0 -r 100 -s | sortBed | gzip -c > ${OUTDIR}/${EXPERIMENT_NAME}.bwa.ext.bed

    fi

}

# assumes fastq format
for NAME in GSM810993_E2F7_Hela GSM810994_Input_Hela
do
    run_bwa_pipeline ${ROOTDIR} ${NAME} ${NAME}.081213 $OUT_DIR 
done

