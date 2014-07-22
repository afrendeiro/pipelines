#!/bin/bash

usage() {
	echo "Usage: $0 -i <input file> -g <genome reference> [-o]"
	echo "-i	Input file [peaks.tsv]"
	echo "-g 	Chromossome sizes file (e.g. )"
	echo "-e	Annotation of genomic regions to exclude in random peaks."
	echo "-l	Genome location annotation file."
	echo "-t	TSS annotation file"
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
while getopts hi:g:e:l:t:o: opt
do
    case "$opt" in
        h)  usage;;
        i)  INFILE="`readlink -m $OPTARG`";;
		g)  GENOME="`readlink -m $OPTARG`";;
		l)  GENOMELOCS="`readlink -m $OPTARG`";;
		e)  EXCLUDE="`readlink -m $OPTARG`";;
		t)	TSSs="`readlink -m $OPTARG`";;
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

BASENAME=`basename $INFILE | sed 's/\(.*\)\..*/\1/'`

BASENAME_SHUFFLED=${BASENAME}.shuffled


# generate shuffled peaks
echo "generating random peaks for $BASENAME"
bedtools shuffle -i $INFILE -g $GENOME -excl $EXCLUDE -seed 0 | bedtools sort > ${OUTDIR}/$BASENAME_SHUFFLED.bed

# get peaks sumits in bed format (not valid for ChIP-chip data)
echo "extracting summits for $BASENAME"	
tail -n +2 $INFILE | awk '{print $1,$5,($5+1)}' | sed '{:q;N;s/ /\t/g;t q}' > ${OUTDIR}/${BASENAME}.summits.bed

echo "extracting summits for $BASENAME_SHUFFLED"	
tail -n +2 ${OUTDIR}/$BASENAME_SHUFFLED.bed | awk '{print $1,$5,($5+1)}' | sed '{:q;N;s/ /\t/g;t q}' > ${OUTDIR}/${BASENAME_SHUFFLED}.summits.bed

# Intersect with genomic locations
echo "overlaping $BASENAME with genomic locations"	
bedtools intersect -wao -a $INFILE -b $GENOMELOCS | cut -f 1,2,3,5,9 > ${OUTDIR}/${BASENAME}.location.bed

echo "overlaping $BASENAME_SHUFFLED with genomic locations"	
bedtools intersect -wao -a ${OUTDIR}/$BASENAME_SHUFFLED.bed -b $GENOMELOCS | cut -f 1,2,3,5,9 > ${OUTDIR}/${BASENAME_SHUFFLED}.location.bed

# Measure distance to TSS
echo "measuring distance to TSS for $BASENAME"	
bedtools closest -d -a $INFILE -b $TSSs | cut -f 1,2,7,9,11 > ${OUTDIR}/${BASENAME}.TSS.bed

echo "measuring distance to TSS for $BASENAME_SHUFFLED"	
bedtools closest -d -a ${OUTDIR}/$BASENAME_SHUFFLED.bed -b $TSSs | cut -f 1,2,7,9,11 > ${OUTDIR}/${BASENAME_SHUFFLED}.TSS.bed



