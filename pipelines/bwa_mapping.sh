#!/bin/bash

if [ $# -lt 3 ] ; then
   echo "You must specify at least 3 arguments. See help with $0 -h"
   exit 1
fi

usage() {
  echo "Usage: $0 -r <path> -i <path> -o <path> -g <file> [-h]"
  echo "-r rootdir"
  exit
}

#Process the arguments
    # options that require arguments have a : in front of them
while getopts hi:o:g: opt
do
    case "$opt" in
        h)  usage;;
        i)  INFILE="`readlink -e $OPTARG`";;
        o)  OUTDIR="`readlink -f $OPTARG`";;
        \?) usage;;
   esac
done

# Capture additional genome directory argument if given, if not set to current dir
shift $(($OPTIND - 1))
if [[ $1 ]]; then
    GENOMEDIR="`readlink -f $1`"
else
    GENOMEDIR="`readlink -f .`"
fi
echo $GENOMEDIR

SAMPLE=${INFILE%.*}
echo $SAMPLE

exit

GENOMEREF="data/human/annotation/hg19/hg19.fa"
GENOMESIZES="data/human/annotation/hg19/hg19_ChrSizes.tsv"

if [ -s ${OUTDIR}/${INFILE}.bwa.bam ]; then
    echo "Results for $INFILE already exist: ${OUTDIR}/${INFILE}.bwa.bam, skipping.."
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
