#!/bin/bash
 
# Set default values for paths
ROOTDIR="."
ANNOTATIONDIR=${ROOTDIR}"/annotation"
OUTDIR="./output"
 
#Check to see if at least one argument was specified
if [ $# -lt 1 ] ; then
   echo "You must specify at least 1 argument."
   exit 1
fi

usage() {
  echo "Usage: $0 [-r <path>] [-a <path>] [-o <path>] [-h]"
}

#Process the arguments
    # options that require arguments have a : in front of them
while getopts hr:a:o: opt
do
   case "$opt" in
      h) usage;;
      r) ROOTDIR=$OPTARG;;
      a) ANNOTATIONDIR=$OPTARG;;
      o) OUTDIR=$OPTARG;;
      \?) usage;;
   esac
done

echo $ROOTDIR
echo $ANNOTATIONDIR
echo $OUTDIR

