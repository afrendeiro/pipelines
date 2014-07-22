
PEAKS='output/oikopeura/ChIP-chip/ChIP-chip_Male_E2F7.peaks.bed'
GENOME='data/oikopleura/assembly/Oikopleura_reference_chrSizes.tsv'
REPEATSGAPS='data/oikopleura/annotation/repeats_gaps_merged.bed'
OUTPUT='output/oikopeura/ChIP-chip'

SAMPLENAME=`basename $PEAKS | sed -e 's/\..*$//'`

bedtools shuffle -i $PEAKS -g $GENOME -excl $REPEATSGAPS -seed 0 | bedtools sort > ${OUTPUT}/${SAMPLENAME}.peaks.shuffled.bed
