#!/bin/bash

usage() {
	echo "Usage: $0 -i <input file> -g <genome windows directory> [-o]"
	echo "-i	Input file in bed format"
	echo "-w	Directory with genomic windows (e.g. data/human/genome_windows/TSS)."
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
while getopts hi:w:o: opt
do
    case "$opt" in
        h)  usage;;
        i)  INFILE="`readlink -m $OPTARG`";;
		w)  WINDOWSDIR="`readlink -m $OPTARG`";;
		o)	OUTDIR="`readlink -m $OPTARG`"
			mkdir -p $OUTDIR;;
        \?) usage;;
   esac
done

# If output directory argument is not set, set to current dir
shift $(($OPTIND - 1))
if [[ ! $OUTDIR ]]; then
	OUTDIR="`readlink -f ./`"
fi

BASENAME=`basename $INFILE .bwa.bed`

WINDOWS=`ls $WINDOWSDIR/heatmap_models_*`

# Function to calculate coverage in those windows for each ChIP sample
for WINDOW in $WINDOWS
do
	WINDOWNAME=${WINDOW/heatmap_models_/}
	WINDOWNAME=`basename ${WINDOWNAME/.txt/}`
	echo "doing $BASENAME sample on window $WINDOWNAME"
	intersectBed -wa -c -a $WINDOW -b $INFILE | sortBed | cut -f 1,2,3,4,8 > $OUTDIR/${BASENAME}.$WINDOWNAME.bed
done
