#! /bin/bash

# Function to calculate coverage in those windows for each ChIP sample
calculate_coverage () {
	WINDOWSDIR=$2
	WINDOWS=`ls $WINDOWSDIR`
	REFDIR=$3
	OUTDIR=$4
	sample=$1
	SAMPLE=${REFDIR}/$sample
	for window in $WINDOWS
	do
	echo "doing ${sample/_*} sample on window $WINDOW"
	WINDOW=${window/heatmap_models_/}
	WINDOW=${WINDOW/.txt/}
	
	intersectBed -wa -wb -a $SAMPLE -b ${WINDOWSDIR}/$window | sortBed | cut -f 1,2,3,4,8 > $OUTDIR/${sample/_*}.$WINDOW.TTScoverage.bed
	done
}

# Run on samples for coverage:
# add full path to sample file (extended reads bed file)
REFDIR="/home/s3/afr/data/human/E2F7"
OUTDIR="/home/s3/afr/data/human/E2F7/TSS_coverage"
SAMPLES="E2F7_genome_coverage_norm.bed"

## for TSS coverage

for sample in $SAMPLES
do
	calculate_coverage $sample ~/data/human/genome_windows/TSS/ $REFDIR $OUTDIR
done

REFDIR="/home/s3/afr/data/human/E2F7"
OUTDIR="/home/s3/afr/data/human/E2F7/TTS_coverage"
SAMPLES="E2F7_genome_coverage_norm.bed"
## for TTS coverage
for sample in $SAMPLES
do
	calculate_coverage $sample ~/data/human/genome_windows/TTS/ $REFDIR $OUTDIR
done