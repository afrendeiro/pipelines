#!/bin/bash

usage() {
	echo "Usage: $0 -i <input file> -g <genome reference> [-o]"
	echo "-i	Input file [peaks.tsv]"
	echo "-g	Genome reference name. Should not have extension (e.g. data/human/hg19/hg19)."
	echo "-d	MEME database file to use."
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
while getopts hi:g:d:o: opt
do
    case "$opt" in
        h)  usage;;
        i)  INFILE="`readlink -m $OPTARG`";;
		g)  GENOMEREF=$OPTARG;;
		d)	MEMEDB=$OPTARG;;
		o)	OUTDIR="`readlink -m $OPTARG`"
			mkdir -p $OUTDIR;;
        \?) usage;;
   esac
done

# If output directory argument is not set, set to current dir
shift $(($OPTIND - 1))
if [[ ! $OUTDIR ]]; then
	OUTDIR="`readlink -m ./`"
fi

BASENAME=`basename $INFILE .peaks.bed`

# get peaks sumits in bed format
tail -n +2 $INFILE | awk '{print $1,$5,($5+1)}' | sed '{:q;N;s/ /\t/g;t q}' > ${OUTDIR}/${BASENAME}.peaks.summits.bed

# get 500bp region around peak summit and get fasta sequences
bedtools slop -i ${OUTDIR}/${BASENAME}.peaks.summits.bed -g ${GENOMEREF}_chrsizes.tsv -b 250 | bedtools getfasta -fi ${GENOMEREF}.fa -bed - -fo ${OUTDIR}/{$BASENAME}.peaks.fa

# small curated database
meme-chip -meme-p 24 -db $MEMEDB -oc ${OUTDIR}/meme-chip ${OUTDIR}/{$BASENAME}.peaks.fa