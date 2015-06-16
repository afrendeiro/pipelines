#!/usr/bin/env python


def slurmHeader(jobName, output, queue="shortq", ntasks=1, time="10:00:00",
                cpusPerTask=16, memPerCpu=2000, nodes=1, userMail=""):
    cmd = """    #!/bin/bash
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

    return cmd


def slurmFooter():
    cmd = "    date"

    return cmd


def slurmSubmitJob(jobFile):
    import os
    cmd = "sbatch %s" % jobFile

    return os.system(cmd)


def removeFile(fileName):
    cmd = "rm {0}".format(fileName)

    return cmd


def moveFile(old, new):
    cmd = "mv {0} {1}".format(old, new)

    return cmd


def makeDir(directory):
    cmd = "mkdir -p {0}".format(directory)

    return cmd


def mergeBams(inputBams, outputBam):
    cmd = "java -Xmx4g -jar /cm/shared/apps/picard-tools/1.118/MergeSamFiles.jar"
    cmd += " USE_THREADING=TRUE"
    cmd += " " + (" ".join(["INPUT=%s"] * len(inputBams))) % tuple(inputBams)
    cmd += " OUTPUT={0}".format(outputBam)

    return cmd


def fastqc(inputBam, outputDir, sampleName):
    import os

    initial = os.path.splitext(os.path.basename(inputBam))[0]

    cmd1 = "fastqc --noextract --outdir {0} {1}".format(outputDir, inputBam)

    cmd2 = "mv {0}_fastqc.html {1}_fastqc.html".format(os.path.join(outputDir, initial), os.path.join(outputDir, sampleName))

    cmd3 = "mv {0}_fastqc.zip {1}_fastqc.zip".format(os.path.join(outputDir, initial), os.path.join(outputDir, sampleName))

    return [cmd1, cmd2, cmd3]


def bam2fastq(inputBam, outputFastq, outputFastq2=None, unpairedFastq=None):
    cmd = "java -Xmx4g -jar /cm/shared/apps/picard-tools/1.118/SamToFastq.jar"
    cmd += " INPUT={0}".format(inputBam)
    cmd += " FASTQ={0}".format(outputFastq)
    if outputFastq2 is not None and unpairedFastq is not None:
        cmd += " SECOND_END_FASTQ={0}".format(outputFastq2)
        cmd += " UNPAIRED_FASTQ={0}".format(unpairedFastq)

    return cmd


def trimmomatic(inputFastq1, outputFastq1, cpus, adapters, log,
                inputFastq2=None, outputFastq1unpaired=None,
                outputFastq2=None, outputFastq2unpaired=None):

    PE = False if inputFastq2 is None else True
    pe = "PE" if PE else "SE"

    cmd = "java -Xmx4g -jar `which trimmomatic-0.32.jar`"
    cmd += " {0} -threads {1} -trimlog {2} {3}".format(pe, cpus, log, inputFastq1)
    if PE:
        cmd += " {0}".format(inputFastq2)
    cmd += " {0}".format(outputFastq1)
    if PE:
        cmd += " {0} {1} {2}".format(outputFastq1unpaired, outputFastq2, outputFastq2unpaired)
    cmd += " ILLUMINACLIP:{0}:1:40:15:8:true".format(adapters)
    cmd += " LEADING:3 TRAILING:3"
    cmd += " SLIDINGWINDOW:4:10"
    cmd += " MINLEN:36"

    return cmd


def skewer(inputFastq1, outputPrefix, cpus, adapters, inputFastq2=None):

    PE = False if inputFastq2 is None else True
    mode = "pe" if PE else "any"

    cmd = "skewer --quiet"
    cmd += " -t {0}".format(cpus)
    cmd += " -m {0}".format(mode)
    cmd += " -x {0}".format(adapters)
    cmd += " -o {0}".format(outputPrefix)
    cmd += " {0}".format(inputFastq1)
    if inputFastq2 is not None:
        cmd += " {0}".format(inputFastq2)

    return cmd


def bowtie2Map(inputFastq1, outputBam, log, metrics, genomeIndex, maxInsert, cpus, inputFastq2=None):
    import re

    outputBam = re.sub("\.bam$", "", outputBam)
    # Admits 2000bp-long fragments (--maxins option)
    cmd = "bowtie2 --very-sensitive -p {0}".format(cpus)
    cmd += " -x {0}".format(genomeIndex)
    cmd += " --met-file {0}".format(metrics)
    if inputFastq2 is None:
        cmd += " {0} ".format(inputFastq1)
    else:
        cmd += " --maxins {0}".format(maxInsert)
        cmd += " -1 {0}".format(inputFastq1)
        cmd += " -2 {0}".format(inputFastq2)
    cmd += " 2> {0} | samtools view -S -b - | samtools sort - {1}".format(log, outputBam)

    return cmd


def topHatMap(inputFastq, outDir, genome, transcriptome, cpus):
    # TODO:
    # Allow paired input
    cmd = "tophat --GTF {0} --b2-L 15 --library-type fr-unstranded --mate-inner-dist 120".format(transcriptome)
    cmd += " --max-multihits 100 --no-coverage-search"
    cmd += " --num-threads {0} --output-dir {1} {2} {3}".format(cpus, outDir, genome, inputFastq)

    return cmd


def markDuplicates(inputBam, outputBam, metricsFile, tempDir="."):
    import re

    transientFile = re.sub("\.bam$", "", outputBam) + ".dups.nosort.bam"
    outputBam = re.sub("\.bam$", "", outputBam)

    cmd1 = "java -Xmx4g -jar `which MarkDuplicates.jar`"
    cmd1 += " INPUT={0}".format(inputBam)
    cmd1 += " OUTPUT={0}".format(transientFile)
    cmd1 += " METRICS_FILE={0}".format(metricsFile)
    cmd1 += " VALIDATION_STRINGENCY=LENIENT"
    cmd1 += " TMP_DIR={0}".format(tempDir)

    # Sort bam file with marked duplicates
    cmd2 = "samtools sort {0} {1}".format(transientFile, outputBam)

    # Remove transient file
    cmd3 = "if [[ -s {0} ]]; then rm {0}; fi".format(transientFile)

    return [cmd1, cmd2, cmd3]


def removeDuplicates(inputBam, outputBam, cpus=16):
    cmd = "sambamba markdup -t {0} -r {1} {2}".format(cpus, inputBam, outputBam)

    return cmd


def shiftReads(inputBam, genome, outputBam):
    import re

    outputBam = re.sub("\.bam$", "", outputBam)

    cmd = "samtools view -h {0} |".format(inputBam)
    cmd += " shift_reads.py {0} |".format(genome)
    cmd += " samtools view -S -b - |"
    cmd += " samtools sort - {0}".format(outputBam)

    return cmd


def sortIndexBam(inputBam, outputBam):
    import re

    tmpBam = re.sub("\.bam", ".sorted", inputBam)

    cmd1 = "samtools sort {0} {1}".format(inputBam, tmpBam)

    cmd2 = "mv {0}.bam {1}".format(tmpBam, outputBam)

    cmd3 = "samtools index {0}".format(outputBam)

    return [cmd1, cmd2, cmd3]


def indexBam(inputBam):
    cmd = "samtools index {0}".format(inputBam)

    return cmd


def peakTools(inputBam, output, plot, cpus):
    cmd = "Rscript `which run_spp.R` -rf -savp -savp={0} -s=0:5:500 -c={1} -out={2}".format(plot, inputBam, output)

    return cmd


def qc():
    """
    $PRESEQ c_curve -B $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
    -o $PROJECTDIR/mapped/$SAMPLE_NAME_qc_c_curve.txt

    $PRESEQ lc_extrap -e 1e8 -s 2e6 -B $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
    -o $PROJECTDIR/mapped/$SAMPLE_NAME_qc_lc_extrap.txt
    """
    raise NotImplementedError("Function not implemented yet.")


def bamToBigWig(inputBam, outputBigWig, genomeSizes, genome, tagmented=False):
    import os
    import re

    # TODO:
    # addjust fragment length dependent on read size and real fragment size
    # (right now it asssumes 50bp reads with 180bp fragments)

    transientFile = os.path.abspath(re.sub("\.bigWig", "", outputBigWig))

    cmd1 = "bedtools bamtobed -i {0} |".format(inputBam)
    if not tagmented:
        cmd1 += " bedtools slop -i stdin -g {0} -s -l 0 -r 130 |".format(genomeSizes)
        cmd1 += " fix_bedfile_genome_boundaries.py {0} |".format(genome)
    else:
        cmd1 = "bedtools bamtobed -i {0} |".format(inputBam)
        cmd1 += " get5primePosition.py |"
    cmd1 += " genomeCoverageBed -i stdin -bg -g {0} > {1}.cov".format(genomeSizes, transientFile)

    cmd2 = "bedGraphToBigWig {0}.cov {1} {2}".format(transientFile, genomeSizes, outputBigWig)
    # remove cov file
    cmd3 = "if [[ -s {0}.cov ]]; then rm {0}.cov; fi".format(transientFile)

    cmd4 = "chmod 755 {0}".format(outputBigWig)

    return [cmd1, cmd2, cmd3, cmd4]


def addTrackToHub(sampleName, trackURL, trackHub, colour, fivePrime=""):
    cmd1 = """echo "track type=bigWig name='{0} {1}' description='{0} {1}'""".format(sampleName, fivePrime)
    cmd1 += """ height=32 visibility=full maxHeightPixels=32:32:25 bigDataUrl={0} color={1}" >> {2}""".format(trackURL, colour, trackHub)

    cmd2 = "chmod 755 {0}".format(trackHub)

    return [cmd1, cmd2]


def linkToTrackHub(trackHubURL, fileName, genome):
    import textwrap

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


def htSeqCount(inputBam, gtf, output):
    sam = inputBam.replace("bam", "sam")

    cmd1 = "samtools view {0} > {1}".format(inputBam, sam)

    cmd2 = "htseq-count -f sam -t exon -i transcript_id -m union {0} {1} > {2}".format(sam, gtf, output)

    cmd3 = "rm {0}".format(sam)

    return [cmd1, cmd2, cmd3]


def kallisto(inputFastq1, outputDir, transcriptomeIndex, inputFastq2=None):
    cmd = "kallisto quant --plaintext  -l 180 -i {0} -o {1} {2}".format(transcriptomeIndex, outputDir, inputFastq1)
    if inputFastq2 is not None:
        cmd += " {0}".format(inputFastq2)
    return cmd


def genomeWideCoverage(inputBam, genomeWindows, output):
    cmd = "bedtools coverage -abam -counts -a {0} -b {1} > {2}".format(inputBam, genomeWindows, output)

    return cmd


def calculateFRiP(inputBam, inputBed, output):
    cmd = "cut -f 1,2,3 {0} |".format(inputBed)
    cmd += " bedtools coverage -counts -abam {0} -b - |".format(inputBam)
    cmd += " awk '{{sum+=$4}} END {{print sum}}' > {0}".format(output)

    return cmd


def macs2CallPeaks(treatmentBam, outputDir, sampleName, genome, controlBam=None, broad=False):
    if not broad:
        cmd = "macs2 callpeak -t {0}".format(treatmentBam)
        if controlBam is not None:
            cmd += " -c {0}".format(controlBam)
        cmd += " --bw 200 -g {0} -n {0} --outdir {0}".format(genome, sampleName, outputDir)
        # --fix-bimodal --extsize 180
    else:
        # Parameter setting for broad factors according to Nature Protocols (2012)
        # Vol.7 No.9 1728-1740 doi:10.1038/nprot.2012.101 Protocol (D) for H3K36me3
        cmd = "macs2 callpeak -t {0}".format(treatmentBam)
        if controlBam is not None:
            cmd += " -c {0}".format(controlBam)
        cmd += " --broad --nomodel --extsize 73 --pvalue 1e-3 -g {0} -n {1} --outdir {2}".format(
            genome, sampleName, outputDir
        )

    return cmd


def macs2PlotModel(sampleName, outputDir):
    import os

    # run macs r script
    cmd1 = "Rscript {0}/{1}_model.r".format(os.getcwd(), sampleName)
    # move to sample dir
    cmd2 = "mv {0}/{1}_model.pdf {2}/{1}_model.pdf".format(os.getcwd(), sampleName, outputDir)

    return [cmd1, cmd2]


def sppCallPeaks(treatmentBam, controlBam, treatmentName, controlName, outputDir, broad, cpus):
    broad = "TRUE" if broad else "FALSE"

    cmd = "Rscript `which spp_peak_calling.R` {0} {1} {2} {3} {4} {5} {6}""".format(
        treatmentBam, controlBam, treatmentName, controlName, broad, cpus, outputDir
    )

    return cmd


def bamToBed(inputBam, outputBed):
    cmd = "bedtools bamtobed -i {0} > {1}".format(inputBam, outputBed)

    return cmd


def zinbaCallPeaks(treatmentBed, controlBed, cpus, tagmented=False):
    fragmentLength = 80 if tagmented else 180

    cmd = "Rscript `which zinba.R` -l {0} -t {1} -c {2}".format(fragmentLength, treatmentBed, controlBed)

    return cmd


def homerFindMotifs(peakFile, genome, outputDir, size=150, length="8,10,12,14,16", n_motifs=12):

    cmd = "findMotifsGenome.pl {0} {1} {2}".format(peakFile, genome, outputDir)
    cmd += " -mask -size {0} -len {1} -S {2}".format(size, length, n_motifs)

    return cmd


def AnnotatePeaks(peakFile, genome, motifFile, outputBed):
    cmd = "annotatePeaks.pl {0} {1} -mask -mscore -m {2} |".format(peakFile, genome, motifFile)
    cmd += "tail -n +2 | cut -f 1,5,22 > {3}".format(outputBed)

    return cmd


def centerPeaksOnMotifs(peakFile, genome, windowWidth, motifFile, outputBed):

    cmd = "annotatePeaks.pl {0} {1} -size {2} -center {3} |".format(peakFile, genome, windowWidth, motifFile)
    cmd += " awk -v OFS='\t' '{print $2, $3, $4, $1, $6, $5}' |"
    cmd += """ awk -v OFS='\t' -F '\t' '{ gsub("0", "+", $6) ; gsub("1", "-", $6) ; print }' |"""
    cmd += " fix_bedfile_genome_boundaries.py {0} | sortBed > {1}".format(genome, outputBed)

    return cmd


def peakAnalysis(inputBam, peakFile, plotsDir, windowWidth, fragmentsize,
                 genome, n_clusters, strand_specific, duplicates):
    import os

    cmd = "python {0}/lib/peaks_analysis.py {1} {2} {3}".format(
        os.path.abspath(os.path.dirname(os.path.realpath(__file__))),
        inputBam, peakFile, plotsDir
    )
    cmd += " --window-width {4} --fragment-size {5} --genome {6} --n_clusters {7}".format(
        windowWidth, fragmentsize, genome, n_clusters
    )
    if strand_specific:
        cmd += " --strand-specific "
    if duplicates:
        cmd += " --duplicates"

    return cmd


def tssAnalysis(inputBam, tssFile, plotsDir, windowWidth, fragmentsize, genome,
                n_clusters, strand_specific, duplicates):
    import os

    cmd = "python {0}/lib/tss_analysis.py {1} {2} {3}".format(
        os.path.abspath(os.path.dirname(os.path.realpath(__file__))),
        inputBam, tssFile, plotsDir
    )
    cmd += " --window-width {0} --fragment-size {1} --genome {2} --n_clusters {3}".format(
        windowWidth, fragmentsize, genome, n_clusters
    )
    if strand_specific:
        cmd += " --strand-specific"
    if duplicates:
        cmd += " --duplicates"

    return cmd


def footprintAnalysis():
    raise NotImplementedError("Function not implemented yet.")


def plotCorrelations(inputCoverage, plotsDir):
    import os

    cmd = "python {2}/lib/correlations.py {0} {1}".format(
        plotsDir,
        " ".join(["%s"] * len(inputCoverage)) % tuple(inputCoverage),
        os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    )

    return cmd


def diffBind(inputCSV, jobName, plotsDir):
    import os
    cmd = "Rscript {3}/lib/diffBind_analysis.R {0} {1} {2}".format(
        inputCSV, jobName, plotsDir,
        os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    )

    return cmd
