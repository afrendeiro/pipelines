#!/bin/bash

usage() {
	echo "Usage: $0 -c <ChIP file> -i <input file> [-o]"
	echo "-c	ChIP sample file in bed format"
	echo "-i	Input control file in bed format"
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
while getopts hc:i:o: opt
do
    case "$opt" in
        h)  usage;;
		c)  CHIP="`readlink -e $OPTARG`";;
        i)  INPUT="`readlink -e $OPTARG`";;
		o)	OUTDIR="`readlink -f $OPTARG`";;
        \?) usage;;
   esac
done

# If output directory argument is not set, set to current dir
shift $(($OPTIND - 1))
if [[ ! $OUTDIR ]]; then
	OUTDIR="`readlink -f ./`"
fi

BASENAME=`basename $CHIP .bed`

peakzilla.py $CHIP $INPUT > ${OUTDIR}/$BASENAME.peaks.tsv
cut -f 1,2,3 ${OUTDIR}/$BASENAME.peaks.tsv | tail -n +2 > ${OUTDIR}/$BASENAME.peaks.bed

