#!/usr/bin/env python

"""
Ngs toolkit
Python wrapper functions to bioinformatics programs
"""

import os
import re
import textwrap


def slurmHeader(jobName, output, queue="shortq", ntasks=1, time="10:00:00",
                cpusPerTask=16, memPerCpu=2000, nodes=1, userMail=""):
    command = """    #!/bin/bash
    #SBATCH --partition={0}
    #SBATCH --ntasks={1}
    #SBATCH --time={2}

    #SBATCH --cpus-per-task={3}
    #SBATCH --mem-per-cpu={4}
    #SBATCH --nodes={5}

    #SBATCH --job-name={6}
    #SBATCH --output={7}

    #SBATCH --mail-type=end
    #SBATCH --mail-user={8}

    # Start running the job
    hostname
    date

    """.format(queue, ntasks, time, cpusPerTask, memPerCpu,
               nodes, jobName, output, userMail)

    return command


def slurmFooter():
    command = """
    # Job end
    date

    """

    return command


def slurmSubmitJob(jobFile):
    command = "sbatch %s" % jobFile

    return os.system(command)


def removeFile(fileName):
    command = """
    # Removing file
    rm {0}
    """.format(fileName)

    return command


def moveFile(old, new):
    command = """
    # Moving file

    mv {0} {1}
    """.format(old, new)

    return command


def makeDir(directory):
    command = """
    # Removing file

    mkdir -p {0}
    """.format(directory)

    return command


def mergeBams(inputBams, outputBam):
    command = """
    # Merging bam files from replicates
    echo "Merging bam files from replicates"

    java -Xmx4g -jar /cm/shared/apps/picard-tools/1.118/MergeSamFiles.jar \\
    USE_THREADING=TRUE \\
    {1} \\
    OUTPUT={0}
    """.format(outputBam, (" ".join(["INPUT=%s"] * len(inputBams))) % tuple(inputBams))

    return command


def fastqc(inputBam, outputDir):
    command = """
    # Measuring sample quality with Fastqc
    echo "Measuring sample quality with Fastqc"
    module load java/jdk/1.7.0_65
    module load FastQC/0.11.2

    fastqc --noextract \\
    --outdir {0} \\
    {1}

    """.format(outputDir, inputBam)

    return command


def bam2fastq(inputBam, outputFastq, outputFastq2=None, unpairedFastq=None):
    command = """
    # Convert to Fastq format
    echo "Converting to Fastq format"

    java -Xmx4g -jar /cm/shared/apps/picard-tools/1.118/SamToFastq.jar \\
    INPUT={0} \\
    """.format(inputBam)
    if outputFastq2 is None and unpairedFastq is None:
        command += """FASTQ={0}

    """.format(outputFastq)
    else:
        command += """FASTQ={0} \\
    SECOND_END_FASTQ={1} \\
    UNPAIRED_FASTQ={2}

    """.format(outputFastq, outputFastq2, unpairedFastq)

    return command


def trimmomatic(inputFastq1, outputFastq1, cpus, adapters, log,
                inputFastq2=None, outputFastq1unpaired=None,
                outputFastq2=None, outputFastq2unpaired=None):

    PE = False if inputFastq2 is None else True
    pe = "PE" if PE else "SE"

    command = """
    # Trimming adapters from sample
    echo "Trimming adapters from sample"
    module load trimmomatic/0.32

    java -Xmx4g -jar `which trimmomatic-0.32.jar` """

    command += """{0} \\
    -threads {1} \\
    -trimlog {2} \\
    {3} \\
    """.format(pe, cpus, log, inputFastq1)
    if PE:
        command += """\\
    {0} \\
    """.format(inputFastq2)
    command += """\\
    {0} \\
    """.format(outputFastq1)
    if PE:
        command += """\\
    {0} \\
    {1} \\
    {2} \\
    """.format(outputFastq1unpaired,
               outputFastq2, outputFastq2unpaired)
    command += """\\
    ILLUMINACLIP:{0}:1:40:15:8:true \\
    LEADING:3 TRAILING:3 \\
    SLIDINGWINDOW:4:10 \\
    MINLEN:36

    """.format(adapters)

    return command


def skewer(inputFastq1, outputPrefix, cpus, adapters, inputFastq2=None):

    PE = False if inputFastq2 is None else True
    mode = "pe" if PE else "any"

    command = """
    # Trimming adapters from sample
    echo "Trimming adapters from sample"

    skewer \\
    -t {0} \\
    -m {1} \\
    -x {2} \\
    -o {3} \\
    {4} """.format(cpus, mode, adapters, outputPrefix, inputFastq1)

    if inputFastq2 is not None:
        command += """\\
    {0}

    """.format(inputFastq2)

    return command


def bowtie2Map(inputFastq1, outputBam, log, genomeIndex, maxInsert, cpus, inputFastq2=None):
    """

    """
    outputBam = re.sub("\.bam$", "", outputBam)
    # Admits 2000bp-long fragments (--maxins option)
    command = """
    # Mapping reads with Bowtie2
    echo "Mapping reads with Bowtie2"
    module load bowtie/2.2.3
    module load samtools

    bowtie2 --very-sensitive -p {0} \\
    -x {1} \\
    --met-file {2} \\
    """.format(cpus, genomeIndex, log)
    if inputFastq2 is None:
        command += "{0} ".format(inputFastq1)
    else:
        command += """--maxins {0} -1 {1} \\
    -2 {2}""".format(maxInsert, inputFastq1, inputFastq2)
    command += """ | \\
    samtools view -S -b - | \\
    samtools sort - {0}

    """.format(outputBam)

    return command


def shiftReads(inputBam, genome, outputBam):
    outputBam = re.sub("\.bam$", "", outputBam)
    # TODO:
    # Implement read shifting with HTSeq or Cython
    command = """
    # Shifting read of tagmented sample
    echo "Shifting read of tagmented sample"
    module load samtools
    module load python

    samtools view -h {0} | \\
    python {3}/lib/shift_reads.py {1} | \\
    samtools view -S -b - | \\
    samtools sort - {2}

    """.format(inputBam, genome, outputBam,
               os.path.abspath(os.path.dirname(os.path.realpath(__file__))))

    return command


def markDuplicates(inputBam, outputBam, metricsFile, tempDir="."):
    transientFile = re.sub("\.bam$", "", outputBam) + ".dups.nosort.bam"
    outputBam = re.sub("\.bam$", "", outputBam)
    command = """
    # Marking duplicates with piccard
    echo "Mark duplicates with piccard"
    module load samtools

    java -Xmx4g -jar /cm/shared/apps/picard-tools/1.118/MarkDuplicates.jar \\
    INPUT={0} \\
    OUTPUT={1} \\
    METRICS_FILE={2} \\
    VALIDATION_STRINGENCY=LENIENT \\
    TMP_DIR={3}

    # Sort bam file with marked duplicates
    samtools sort {1} \\
    {4}

    if [[ -s {1} ]]
        then
        rm {1}
    fi

    """.format(inputBam, transientFile, metricsFile, tempDir, outputBam)

    return command


def removeDuplicates(inputBam, outputBam, cpus=16):
    command = """
    # Removing duplicates with sambamba
    echo "Remove duplicates with sambamba"
    sambamba markdup -t {2} -r \\
    {0} \\
    {1}

    """.format(inputBam, outputBam, cpus)

    return command


def indexBam(inputBam):
    command = """
    # Indexing bamfile with samtools
    echo "Indexing bamfile with samtools"
    module load samtools

    samtools index {0}

    """.format(inputBam)

    return command


def qc():
    """
    $PRESEQ c_curve -B $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
    -o $PROJECTDIR/mapped/$SAMPLE_NAME_qc_c_curve.txt

    $PRESEQ lc_extrap -e 1e8 -s 2e6 -B $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
    -o $PROJECTDIR/mapped/$SAMPLE_NAME_qc_lc_extrap.txt
    """
    raise NotImplementedError("Function not implemented yet.")


def bamToUCSC(inputBam, outputBigWig, genomeSizes, genome, tagmented=False):
    transientFile = os.path.abspath(re.sub("\.bigWig", "", outputBigWig))
    if not tagmented:
        # TODO:
        # addjust fragment length dependent on read size and real fragment size
        # (right now it asssumes 50bp reads with 180bp fragments)
        command = """
    # Making bigWig tracks from bam file
    echo "making bigWig tracks from bam file"
    module load bedtools

    bedtools bamtobed -i {0} | \\
    bedtools slop -i stdin -g {1} -s -l 0 -r 130 | \\
    python {5}/lib/fix_bedfile_genome_boundaries.py {4} | \\
    genomeCoverageBed -i stdin -bg -g {1} \\
    > {2}.cov

    bedGraphToBigWig {2}.cov \\
    {1} \\
    {3}

    # remove cov file
    if [[ -s {2}.cov ]]
        then
        rm {2}.cov
    fi

    chmod 755 {3}

    """.format(inputBam, genomeSizes, transientFile, outputBigWig, genome,
               os.path.abspath(os.path.dirname(os.path.realpath(__file__))))
    else:
        command = """
    # Making bigWig tracks from bam file
    echo "making bigWig tracks from bam file"
    module load bedtools

    bedtools bamtobed -i {0} | \\
    python {4}/lib/get5primePosition.py | \\
    genomeCoverageBed -i stdin -bg -g {1} \\
    > {2}.cov

    bedGraphToBigWig {2}.cov \\
    {1} \\
    {3}

    # remove cov file
    if [[ -s {2}.cov ]]
        then
        rm {2}.cov
    fi

    chmod 755 {3}

    """.format(inputBam, genomeSizes, transientFile, outputBigWig,
               os.path.abspath(os.path.dirname(os.path.realpath(__file__))))

    return command


def addTrackToHub(sampleName, trackURL, trackHub, colour, fivePrime=""):
    command = """
    # Adding track to TrackHub
    echo "Adding track to TrackHub"
    echo "track type=bigWig name='{0} {1}' description='{0} {1}' """.format(sampleName, fivePrime)
    command += """height=32 visibility=full maxHeightPixels=32:32:25 bigDataUrl={0} color={1}" >> \\
    {2}

    chmod 755 {2}

    """.format(trackURL, colour, trackHub)

    return command


def linkToTrackHub(trackHubURL, fileName, genome):
    db = "org" if genome == "hg19" else "db"  # different database call for human
    genome = "human" if genome == "hg19" else genome  # change hg19 to human

    html = """
    <html>
        <head>
            <meta http-equiv="refresh" content="0; url=http://genome.ucsc.edu/cgi-bin/hgTracks?"""
    html += """{db}={genome}&hgt.customText={trackHubURL}" />
        </head>
    </html>
    """.format(trackHubURL=trackHubURL, genome=genome, db=db)

    with open(fileName, 'w') as handle:
        handle.write(textwrap.dedent(html))


def genomeWideCoverage(inputBam, genomeWindows, output):
    command = """
    # Calculate genome-wide coverage
    echo "Calculating genome-wide coverage"
    module load bedtools

    bedtools coverage -abam -counts \\
    -a {0} \\
    -b {1} \\
    > {2}

    """.format(inputBam, genomeWindows, output)

    return command


def macs2CallPeaks(treatmentBam, controlBam, outputDir, sampleName, genome, broad=False, plot=True):
    if genome == "hg19":
        genome = "hs"
    elif genome == "mm10":
        genome = "mm"
    elif genome == "dr7":
        raise ValueError("Genome dr7 not yet supported for peak calling with MACS2.")

    if not broad:
        command = """
    # Call peaks with MACS2
    echo "Calling peaks with MACS2"

    macs2 callpeak \\
    -t {0} \\
    -c {1} \\
    --bw 200 \\
    -g {2} -n {3} --outdir {4}

    """.format(treatmentBam, controlBam, genome, sampleName, outputDir)
        command = """
    # Call peaks with MACS2
    echo "Calling peaks with MACS2"

    macs2 callpeak \\
    -t {0} \\
    -c {1} \\
    --bw 200 \\
    --fix-bimodal --extsize 200 \\
    -g {2} -n {3} --outdir {4}

    """.format(treatmentBam, controlBam, genome, sampleName, outputDir)

    else:
        # Parameter setting for broad factors according to Nature Protocols (2012)
        # Vol.7 No.9 1728-1740 doi:10.1038/nprot.2012.101 Protocol (D) for H3K36me3
        command = """
    # Call peaks with MACS2
    echo "Calling peaks with MACS2"

    macs2 callpeak \\
    -t {0} \\
    -c {1} \\
    --broad --nomodel --extsize 73 --pvalue 1e-3 \\
    -g {2} -n {3} --outdir {4}

    """.format(treatmentBam, controlBam, genome, sampleName, outputDir)

    if plot:
        command += """
    # Plot MACS2 crosscorrelation

    echo "Plotting MACS2 crosscorrelation"
    module load R

    Rscript {0}/{1}_model.r

    mv {2}/{1}_model.pdf {0}/{1}_model.pdf

    """.format(outputDir, sampleName, os.getcwd())

    return command


def sppCallPeaks(treatmentBam, controlBam, treatmentName, controlName, outputDir, broad, cpus):
    broad = "TRUE" if broad else "FALSE"
    command = """
    # Calling peaks with SPP
    echo "calling peaks with SPP"
    module load R

    Rscript {6}/lib/spp_peak_calling.R \\
    {0} \\
    {1} \\
    {2} \\
    {3} \\
    {4} \\
    {5} \\
    {6}

    """.format(treatmentBam, controlBam, treatmentName, controlName, broad, cpus,
               os.path.abspath(os.path.dirname(os.path.realpath(__file__))))

    return command


def homerFindMotifs(peakFile, genome, outputDir, size=150, length="8,10,12,14,16", n_motifs=12):
    command = """
    # Find motifs with Homer
    echo "finding motifs with Homer"

    findMotifsGenome.pl {0} \\
    {1} \\
    {2} \\
    -mask -size {3} -len {4} -S {5}

    """.format(peakFile, genome, outputDir, size, length, n_motifs)

    return command


def AnnotatePeaks(peakFile, genome, motifFile, outputBed):
    command = """
    # Annotate peaks with motif score
    echo "annotating peaks with motif score"

    annotatePeaks.pl {0} \\
    {1} \\
    -mask -mscore -m {2} | \\
    tail -n +2 | cut -f 1,5,22 \\
    > {3}
    """.format(peakFile, genome, motifFile, outputBed)

    return command


def centerPeaksOnMotifs(peakFile, genome, windowWidth, motifFile, outputBed):
    command = """
    # Center peaks on motif
    echo "centering peaks on motif"

    annotatePeaks.pl {0} \\
    {1} \\
    -size {2} -center {3} | \\
    awk -v OFS='\\t' '{{print $2, $3, $4, $1, $6, $5}}' | \\
    awk -v OFS='\\t' -F '\t' '{{ gsub("0", "+", $6) ; gsub("1", "-", $6) ; print }}' | \\
    python {5}/lib/fix_bedfile_genome_boundaries.py {1} | \\
    sortBed > {4}
    """.format(peakFile, genome, windowWidth, motifFile, outputBed,
               os.path.abspath(os.path.dirname(os.path.realpath(__file__))))

    return command


def peakAnalysis(inputBam, peakFile, plotsDir, windowWidth, fragmentsize,
                 genome, n_clusters, strand_specific, duplicates):
    command = """
    # Analyse peak profiles
    echo "Analysing peak profiles"

    python {0}/lib/peaks_analysis.py {1} \\
    {2} \\
    {3} \\
    --window-width {4} \\
    --fragment-size {5} \\
    --genome {6} \\
    --n_clusters {7} """.format(
        os.path.abspath(os.path.dirname(os.path.realpath(__file__))),
        inputBam, peakFile, plotsDir, windowWidth, fragmentsize, genome, n_clusters
    )
    if strand_specific:
        command += "--strand-specific "
    if duplicates:
        command += "--duplicates\n"

    return command


def tssAnalysis(inputBam, tssFile, plotsDir, windowWidth, fragmentsize, genome,
                n_clusters, strand_specific, duplicates):
    command = """
    # Analyse signal over TSSs
    echo "Analysing signal over TSSs"

    python {0}/lib/tss_analysis.py {1} \\
    {2} \\
    {3} \\
    --window-width {4} \\
    --fragment-size {5} \\
    --genome {6} \\
    --n_clusters {7} """.format(
        os.path.abspath(os.path.dirname(os.path.realpath(__file__))),
        inputBam, tssFile, plotsDir, windowWidth, fragmentsize, genome, n_clusters
    )
    if strand_specific:
        command += "--strand-specific "
    if duplicates:
        command += "--duplicates "

    return command


def footprintAnalysis():
    raise NotImplementedError("Function not implemented yet.")


def plotCorrelations(inputCoverage, plotsDir):
    command = """
    # Plot correlations
    echo "plotting correlations"

    python {2}/lib/correlations.py \\
    {0} \\
    {1}

    """.format(plotsDir, 
               " ".join(["%s"] * len(inputCoverage)) % tuple(inputCoverage),
               os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
               )

    return command


def diffBind(inputCSV, jobName, plotsDir):
    command = """
    # Detect differential binding with diffBind
    echo "Detecting differential binding with diffBind"
    module load R

    Rscript {3}/lib/diffBind_analysis.R \\
    {0} \\
    {1} \\
    {2} \\

    """.format(inputCSV, jobName, plotsDir,
               os.path.abspath(os.path.dirname(os.path.realpath(__file__))))

    return command