
ANOTDIR="/home/s3/afr/data/human/annotation/"
REFDIR="/home/s3/afr/data/human/E2F7"
OUTDIR="/home/s3/afr/output/human/E2F7"

# get housekeeping genes
wget -S http://www.tau.ac.il/~elieis/HKG/HK_genes.txt -O - | cut -f 2 > ${ANOTDIR}/HK_genes.txt 