#!/bin/bash

usage() {
	echo "Usage: $0 -i <infile.fastq> -g <genomeref> [-o]"
	echo "-i	Input file in fastq format"
	echo "-g	Genome reference name. Should not have extension (e.g. data/human/hg19/hg19)"
	echo "		bwa index with same name must exist"
	echo "-o	OPTIONAL: Directory for output files. Will be"
	echo "		created if does not exist. Default: current directory"
	exit
}

if [ $# -lt 2 ] ; then
	echo "You must specify at least 2 arguments."
	echo 
	usage
	exit 1
fi

#Process the arguments
    # options that require arguments have a : in front of them
while getopts hi:g: opt
do
    case "$opt" in
        h)  usage;;
        i)  INFILE="`readlink -e $OPTARG`";;
        g)  GENOMEREF=$OPTARG;;
		o)	OUTDIR="`readlink -f $OPTARG`"
			mkdir -p $OUTDIR;;
        \?) usage;;
   esac
done

# If output directory argument is not set, set to current dir
shift $(($OPTIND - 1))
if [[ ! $OUTDIR ]]; then
	OUTDIR="`readlink -f ./`"
fi

BASENAME=`basename $INFILE .fastq`


if [ -s ${OUTDIR}/${BASENAME}.bwa.bam ]; then
    echo "Results for $BASENAME already exist: ${OUTDIR}/${BASENAME}.bwa.bam, skipping.."
else
    echo "Running bwa alignment for $BASENAME.."

    #align to genome, output .bam file
	bwa aln -t 24 $GENOMEREF ${DIR}/${INFILE} | bwa samse $GENOMEREF - ${DIR}/${INFILE} | samtools view -b -q 30 -u -T $GENOMEREF - | samtools sort - ${OUTDIR}/${BASENAME}.bwa
	
    # create index
    echo "Creating index for $BASENAME"
    samtools index ${OUTDIR}/${BASENAME}.bwa.bam
    mv ${OUTDIR}/${BASENAME}.bwa.bam.bai ${OUTDIR}/${BASENAME}.bwa.bai

    # count mapped reads, export to file
    echo "Counting mapped reads for $BASENAME"
    samtools view -c -F 4 ${OUTDIR}/${BASENAME}.bwa.bam | awk -v var=$BASENAME '$2 = $2 FS var' | sed '{:q;N;s/ /\t/g;t q}' > ${OUTDIR}/${BASENAME}.readcount.txt

    # convert to bed
    echo "Converting $BASENAME to bed format"
    bamToBed -i ${OUTDIR}/${BASENAME}.bwa.bam > ${OUTDIR}/${BASENAME}.bwa.bed

    # extend reads
    #echo "Extending reads for $BASENAME"
    #slopBed -g $GENOMESIZES -l 0 -r 100 -s | sortBed | gzip -c > ${OUTDIR}/${BASENAME}.bwa.ext.bed

fi
