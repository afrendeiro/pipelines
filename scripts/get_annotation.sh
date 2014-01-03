#!/bin/bash

if [ $# -lt 1 ] ; then
   echo "You must specify at least 1 argument. See help with $0 -h"
   exit 1
fi

usage() {
  echo "Usage: $0 [-a <path>] -s specie [-h]"
  echo "-r rootdir"
  exit
}

#Process the arguments
    # options that require arguments have a : in front of them
while getopts hs: opt
do
    case "$opt" in
        h)  usage;;
        s)  SPECIE=$OPTARG;;
        \?) usage;;
   esac
done

# Capture additional base directory argument if given, if not set to current dir
echo $3
if [[ $3 ]]; then
    BASEDIR="`readlink -f $3`"
else
    BASEDIR="`readlink -f .`"
fi

ANNOTDIR=${BASEDIR}/data/${SPECIE}/annotation
# Make directory for specie annotation if doesn't exist
mkdir -p $ANNOTDIR

exit

#### GET ANNOTATION based on species ####


# Human genome
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit -O - > ${ANNOTDIR}/hg19.2bit
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/md5sum.txt

# get housekeeping genes
wget -S http://www.tau.ac.il/~elieis/HKG/HK_genes.txt -O - | cut -f 2 > ${ANNOTDIR}/HK_genes.txt 