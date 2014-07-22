#!/bin/bash
# Pipeline for PE samples

# paths and variables to change
RAW=/sysdev/s3/share/data/oikopleura/chip-seq/raw
MAPPED=/sysdev/s3/share/data/oikopleura/chip-seq/mapped

GENOMEREF=~/data/oikopleura/assembly/Oikopleura_reference_unmasked_v3.0.fa
CHRSIZES=~/data/oikopleura/assembly/Oikopleura_reference_chrSizes.tsv

SAMPLES=`echo {'TB','D2','D6-F','D6-M'}_{'IN','E2F1','E2F7','POL2'}_R1`
#SAMPLES=`echo {'TB','D2','D6-F','D6-M'}_{'IN','E2F1','E2F7','POL2'}_{'R1','R2'}`

for SAMPLE in $SAMPLES
do
    # remove adapters, trim reads 
    #if [ ! -s ${RAW}/${SAMPLE}_'S1'.trimmed.fastq.gz ] && [ ! -s ${RAW}/${SAMPLE}_'S2'.trimmed.fastq.gz ]; then
    #    echo `date +%Y.%m.%d-%H:%M:%S` - "Removing adapters and trimming reads for $SAMPLE..." >> log.txt
    #    
    #    # settings optimized for 100bp PE illumina reads
    #    # use alias trimmomatic='java -jar  ~/Trimmomatic-0.32/trimmomatic-0.32.jar'
    #    trimmomatic PE \
    #    ${RAW}/${SAMPLE}_{'S1','S2'}.fastq.gz ${RAW}/${SAMPLE}_{'S1','S2'}.trimmed_{'P','U'}.fastq \
    #    ILLUMINACLIP:/home/s3/afr/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:1:40:15 \ # change path of adaptors if needed
    #    SLIDINGWINDOW:5:15 \
    #    MINLEN:36
    #
    #    # use unpaired reads 
    #    cat ${RAW}/${SAMPLE}_'S1'.trimmed_{'P','U'}.fastq | gzip > ${RAW}/${SAMPLE}_'S1'.trimmed.fastq.gz
    #    cat ${RAW}/${SAMPLE}_'S2'.trimmed_{'P','U'}.fastq | gzip > ${RAW}/${SAMPLE}_'S2'.trimmed.fastq.gz
    #fi

    if [ ! -s ${RAW}/${SAMPLE}.bwa.bam ]; then
        echo `date +%Y.%m.%d-%H:%M:%S` - "Running bwa alignment for $SAMPLE..." >> log.txt
        for SINGLE in ${SAMPLE}_{'S1','S2'}
        do
            # map to genome
            echo `date +%Y.%m.%d-%H:%M:%S` - "Aligning $SINGLE..." >> log.txt
            bwa aln -t 24 $GENOMEREF ${RAW}/$SINGLE.fastq.gz > ${MAPPED}/$SINGLE.sai
        done
    fi

    # produce paired alignment, filter reads Q<30, sort by name and output .bam file
    if [ -s ${MAPPED}/${SAMPLE}_S1.sai ] && [ -s ${MAPPED}/${SAMPLE}_S2.sai ] && [ ! -s ${MAPPED}/${SAMPLE}.bwa.bam ]; then
        echo `date +%Y.%m.%d-%H:%M:%S` - "Producing paired alignment for $SAMPLE..." >> log.txt
        bwa sampe $GENOMEREF ${MAPPED}/${SAMPLE}_{'S1','S2'}.sai ${RAW}/${SAMPLE}_{'S1','S2'}.fastq.gz | samtools view -b -q 30 -u -T $GENOMEREF - | samtools sort -n - ${MAPPED}/${SAMPLE}.bwa
    fi

    # mark PCR duplicates
    if [ -s ${MAPPED}/${SAMPLE}.bwa.bam ]; then
        echo `date +%Y.%m.%d-%H:%M:%S` - "Marking duplicates for $SAMPLE..." >> log.txt
        sambamba markdup -t 12 ${MAPPED}/${SAMPLE}.bwa.bam tmp
        mv tmp ${MAPPED}/${SAMPLE}.bwa.bam
    fi

    # create index
    if [ -s ${MAPPED}/${SAMPLE}.bwa.bam ] && [ ! -s ${MAPPED}/${SAMPLE}.bwa.bam.bai ]; then
        echo `date +%Y.%m.%d-%H:%M:%S` - "Creating index for $SAMPLE..." >> log.txt
        sambamba index -t 12 ${MAPPED}/${SAMPLE}.bwa.bam
    fi

    # count mapped reads, save to file
    if [ -s ${MAPPED}/${SAMPLE}.bwa.bam ]; then
        echo `date +%Y.%m.%d-%H:%M:%S` - "Counting mapped reads for $SAMPLE..." >> log.txt
        sambamba view -c -F "not unmapped" -t 12 ${MAPPED}/${SAMPLE}.bwa.bam | awk -v var=$SAMPLE '$2 = $2 FS var' | sed '{:q;N;s/ /\t/g;t q}' >> ${MAPPED}/readcounts.tsv
    fi

    # compute statistics, save to file
    if [ -s ${MAPPED}/${SAMPLE}.bwa.bam ] && [ ! -s ${MAPPED}/${SAMPLE}.bwa.flagstat.txt ]; then
        echo `date +%Y.%m.%d-%H:%M:%S` - "Computing statistics for $SAMPLE..." >> log.txt
        sambamba flagstat -t 12 ${MAPPED}/${SAMPLE}.bwa.bam > ${MAPPED}/${SAMPLE}.bwa.flagstat.txt
    fi

    # estimate insert size (from only 500k reads)
    if [ -s ${MAPPED}/${SAMPLE}.bwa.bam ] && [ ! -s ${MAPPED}/${SAMPLE}.bwa.insert.txt ]; then
        echo `date +%Y.%m.%d-%H:%M:%S` - "Estimating insert size for $SAMPLE..." >> log.txt
        sambamba view -t 12 ${MAPPED}/${SAMPLE}.bwa.bam | head -n 500000 | getinsertsize.py - > ${MAPPED}/${SAMPLE}.bwa.insert.txt
    fi

    # convert to bed
    if [ -s ${MAPPED}/${SAMPLE}.bwa.bam ] && [ ! -s ${MAPPED}/${SAMPLE}.bwa.bed ]; then
        echo `date +%Y.%m.%d-%H:%M:%S` - "Converting $SAMPLE to bed format..." >> log.txt
        bamToBed -i ${MAPPED}/${SAMPLE}.bwa.bam | sortBed > ${MAPPED}/${SAMPLE}.bwa.bed
        # pe
        #bamToBed -bedpe -i ${MAPPED}/${SAMPLE}.bwa.bam | sortBed > ${MAPPED}/${SAMPLE}.bwa.bed
    fi

    # extend reads, make bedgraph and bigwig 
    if [ -s ${MAPPED}/${SAMPLE}.bwa.bed ] && [ ! -s ${MAPPED}/${SAMPLE}.bwa.bw ]; then
        echo `date +%Y.%m.%d-%H:%M:%S` - "Making Bedgraph and BigWig for $SAMPLE..." >> log.txt
        bedtools slop -s -l 0 -r 99 -i ${MAPPED}/${SAMPLE}.bwa.bed -g $CHRSIZES > ${MAPPED}/${SAMPLE}.bwa.ext.bed
        genomeCoverageBed -bg -i ${MAPPED}/${SAMPLE}.bwa.ext.bed -g $CHRSIZES > ${MAPPED}/${SAMPLE}.bwa.bedgraph
        bedGraphToBigWig ${MAPPED}/${SAMPLE}.bwa.bedgraph $CHRSIZES ${MAPPED}/${SAMPLE}.bwa.bw
    fi
done

cat ${MAPPED}/*insert* > ${MAPPED}/../inserts.txt
cat ${MAPPED}/*flagstat* > ${MAPPED}/../flagstat.txt