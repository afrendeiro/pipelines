#! /bin/bash

# Function to calculate coverage in those windows for each ChIP sample
calculate_coverage () {
	sample=$1
	WINDOWS=$2
	REFDIR=$3
	OUTDIR=$4
	SAMPLE=${REFDIR}/$sample
	echo "doing genome coverage on sample ${sample/.*}"
	coverageBed -a $SAMPLE -b $WINDOWS | sortBed > $OUTDIR/${sample/.*}.genomecoverage.bed
}

# Run on samples for coverage:
# add full path to sample file (extended reads bed file)
REFDIR="/home/s3/afr/data/human/E2F7"
SAMPLES="GSM810993_E2F7_Hela.081213.bwa.ext.bed GSM810994_Input_Hela.081213.bwa.ext.bed"

for sample in $SAMPLES
do
	calculate_coverage $sample ~/data/human/annotation/hg19/200bp_100bpsliding_windows.bed $REFDIR ${REFDIR}/genome_coverage/
done