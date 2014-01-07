#!/bin/bash

usage() {
	echo "Usage: $0 -c <ChIP file> -i <input file> [<outputdir>]"
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
while getopts hi:o:g: opt
do
    case "$opt" in
        h)  usage;;
		c)  CHIP="`readlink -e $OPTARG`";;
        i)  INPUT="`readlink -e $OPTARG`";;
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

BASENAME=`basename $CHIP .bed`


peakzilla.py $CHIP $INPUT > ${OUTDIR}/$BASENAME.peaks.tsv
cut -f 1,2,3 ${OUTDIR}/$BASENAME.peaks.tsv | tail -n +2 > ${OUTDIR}/$BASENAME.peaks.bed

