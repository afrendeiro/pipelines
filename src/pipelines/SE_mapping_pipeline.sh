GENOMEREF=~/data/oikopleura/assembly/Oikopleura_reference_masked_v3.0.fa
for INFILE in {'TB','D2','D6-F','D6-M'}_{'IN','E2F1','E2F7','POL2'}_{'R1','R2'}.fastq.gz
do
    BASENAME=`basename $INFILE .fastq.gz`

    if [ -s ${BASENAME}.bwa.bam ]; then
        echo "Results for $BASENAME already exist: ${BASENAME}.bwa.bam, skipping.."
    else
        echo "Running bwa alignment for $BASENAME..."

        #align to genome, output .bam file
        bwa aln -t 24 $GENOMEREF $INFILE | bwa samse $GENOMEREF - $INFILE | samtools view -b -q 30 -u -T $GENOMEREF - | samtools sort - ${BASENAME}.bwa
        
        # create index
        echo "Creating index for $BASENAME..."
        samtools index ${BASENAME}.bwa.bam

        # count mapped reads, save to file
        echo "Counting mapped reads for $BASENAME..."
        samtools view -c -F 4 ${BASENAME}.bwa.bam | awk -v var=$BASENAME '$2 = $2 FS var' | sed '{:q;N;s/ /\t/g;t q}' >> readcounts.tsv

        # estimate insert size (from only 1M reads)
        samtools view ${BASENAME}.bwa.bam | head -n 1000000 | getinsertsize.py - > ${BASENAME}.bwa.insert.txt

        # convert to bed
        echo "Converting $BASENAME to bed format..."
        bamToBed -i ${BASENAME}.bwa.bam | sortBed > ${BASENAME}.bwa.bed

        # extend reads
        #echo "Extending reads for $BASENAME"
        #slopBed -g $GENOMESIZES -l 0 -r 100 -s | sortBed | gzip -c > ${BASENAME}.bwa.ext.bed

    fi
done