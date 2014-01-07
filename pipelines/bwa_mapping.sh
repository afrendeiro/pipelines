#!/bin/bash

usage() {
	echo "Usage: $0 -i <infile.fastq> -g <genomeref> [<outputdir>]"
	echo "-i	Input file in fastq format"
	echo "-g	Genome reference name. Should not have extension."
	echo "		bwa index with same name must exist"
	echo "-o	OPTIONAL: Directory for output files. Will be"
	echo "		created if does not exist. Default: current directory"
	exit
}

if [ $# -lt 3 ] ; then
	echo "You must specify at least 3 arguments."
	echo 
	usage
	exit 1
fi

#Process the arguments
    # options that require arguments have a : in front of them
while getopts hi:o:g: opt
do
    case "$opt" in
        h)  usage;;
        i)  INFILE="`readlink -e $OPTARG`";;
        g)  GENOMEREF=$OPTARG;;
        \?) usage;;
   esac
done

# Capture additional output directory argument if given, if not set to current dir
shift $(($OPTIND - 1))
if [[ $1 ]]; then
    OUTDIR="`readlink -f $1`"
    mkdir -p $OUTDIR
else
    OUTDIR="`readlink -f .`"
fi

BASENAME=`basename $INFILE .fastq`


if [ -s ${OUTDIR}/${BASENAME}.bwa.bam ]; then
    echo "Results for $BASENAME already exist: ${OUTDIR}/${BASENAME}.bwa.bam, skipping.."
else
    echo "Running bwa alignment for $BASENAME.."

    #align to genome, output .bam file
    # filters for quality >=q30
	bwa aln -t 24 $GENOMEREF ${DIR}/${INFILE} | bwa samse $GENOMEREF - ${DIR}/${INFILE} | samtools view -b -q 30 -u -T $GENOMEREF - | samtools sort - ${OUTDIR}/${BASENAME}.bwa
	
    # create index
    echo "Creating index for $BASENAME"
    samtools index ${OUTDIR}/${BASENAME}.bwa.bam
    mv ${OUTDIR}/${BASENAME}.bwa.bam.bai ${OUTDIR}/${BASENAME}.bwa.bai

    # count mapped reads, export to file
    echo "Counting mapped reads for $BASENAME"
    samtools view -c -F 4 ${OUTDIR}/${BASENAME}.bwa.bam | awk -v var=$BASENAME '$2 = $2 FS var' > ${OUTDIR}/${BASENAME}.readcount.txt

    # convert to bed
    echo "Converting $BASENAME to bed format"
    bamToBed -i ${OUTDIR}/${BASENAME}.bwa.bam > ${OUTDIR}/${BASENAME}.bwa.bed

    # extend reads
    #echo "Extending reads for $BASENAME"
    #slopBed -g $GENOMESIZES -l 0 -r 100 -s | sortBed | gzip -c > ${OUTDIR}/${BASENAME}.bwa.ext.bed

fi
