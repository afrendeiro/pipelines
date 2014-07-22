usage() {
	echo "Usage: $0 -i <input file> -g <genome reference> [-o]"
	echo "-i	Input file [peaks.tsv]"
	echo "-g	Base directory of genome annotation (e.g. data/human/hg19/hg19)"
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
        i)  INFILE="`readlink -m $OPTARG`";;
		g)  GENOMEDIR="`readlink -m $OPTARG`";;
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

BASENAME=`basename $INFILE .peaks.tsv`

exit

# Finds target genes based on TFBS based on proximity to nearest TSS

REFDIR="/home/s3/afr/data/human/E2F7"
OUTDIR="/home/s3/afr/output/human/E2F7"

# Extract peak summit position 
# Find closest gene
# measure distance to TSS
# filter where distance is < X
tail -n +2 ${OUTDIR}/$INFILE | awk '{print $1,$5,($5+1)}' | sed '{:q;N;s/ /\t/g;t q}' | \
bedtools closest -d -a - -b ${GENOMEDIR}/genomic_location/TSS.bed | \
cut -f 1,2,3,7,10 > ${OUTDIR}/E2F7_Hela.targets.bed

cut -f 4 ${OUTDIR}/E2F7_Hela.targets.bed | sed 's/\([^_]\_[^_]*\).*/\1/' > ${OUTDIR}/E2F7_Hela.targets.txt



# Find out if gene is part of Operon
