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
    cmd = "date"

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


def fastqc(inputBam, outputDir):
    cmd1 = "module load java/jdk/1.7.0_65"

    cmd2 = "module load FastQC/0.11.2"

    cmd3 = "fastqc --noextract --outdir {0} {1}".format(outputDir, inputBam)

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

    cmd1 = "module load trimmomatic/0.32"

    cmd2 = "java -Xmx4g -jar `which trimmomatic-0.32.jar`"
    cmd2 += " {0} -threads {1} -trimlog {2} {3}".format(pe, cpus, log, inputFastq1)
    if PE:
        cmd2 += " {0}".format(inputFastq2)
    cmd2 += " {0}".format(outputFastq1)
    if PE:
        cmd2 += " {0} {1} {2}".format(outputFastq1unpaired, outputFastq2, outputFastq2unpaired)
    cmd2 += " ILLUMINACLIP:{0}:1:40:15:8:true".format(adapters)
    cmd2 += " LEADING:3 TRAILING:3"
    cmd2 += " SLIDINGWINDOW:4:10"
    cmd2 += " MINLEN:36"

    return [cmd1, cmd2]


def skewer(inputFastq1, outputPrefix, cpus, adapters, inputFastq2=None):

    PE = False if inputFastq2 is None else True
    mode = "pe" if PE else "any"

    cmd = "skewer --quiet"
    cmd += " -t {0}".format(cpus)
    cmd += " -m {1}".format(mode)
    cmd += " -x {2}".format(adapters)
    cmd += " -o {3}".format(outputPrefix)
    cmd += " {4}".format(inputFastq1)
    if inputFastq2 is not None:
        cmd += " {0}".format(inputFastq2)

    return cmd


def bowtie2Map(inputFastq1, outputBam, log, metrics, genomeIndex, maxInsert, cpus, inputFastq2=None):
    import re

    outputBam = re.sub("\.bam$", "", outputBam)
    # Admits 2000bp-long fragments (--maxins option)
    cmd1 = "module load bowtie/2.2.3"

    cmd2 = "module load samtools"

    cmd3 = "bowtie2 --very-sensitive -p {0}".format(cpus)
    cmd3 += " -x {0}".format(genomeIndex)
    cmd3 += " --met-file {0}".format(metrics)
    if inputFastq2 is None:
        cmd3 += " {0} ".format(inputFastq1)
    else:
        cmd3 += " --maxins {0}".format(maxInsert)
        cmd3 += " -1 {1}".format(inputFastq1)
        cmd3 += " -2 {2}".format(inputFastq2)
    cmd3 += " 2> {0} | samtools view -S -b - | samtools sort - {1}".format(log, outputBam)

    return [cmd1, cmd2, cmd3]


def markDuplicates(inputBam, outputBam, metricsFile, tempDir="."):
    import re

    transientFile = re.sub("\.bam$", "", outputBam) + ".dups.nosort.bam"
    outputBam = re.sub("\.bam$", "", outputBam)

    cmd1 = "module load samtools"

    cmd2 = "java -Xmx4g -jar /cm/shared/apps/picard-tools/1.118/MarkDuplicates.jar"
    cmd2 += " INPUT={0}".format(inputBam)
    cmd2 += " OUTPUT={0}".format(transientFile)
    cmd2 += " METRICS_FILE={0}".format(metricsFile)
    cmd2 += " VALIDATION_STRINGENCY=LENIENT"
    cmd2 += " TMP_DIR={0}".format(tempDir)

    # Sort bam file with marked duplicates
    cmd3 = "samtools sort {1} {4}".format(transientFile, outputBam)

    # Remove transient file
    cmd4 = "if [[ -s {1} ]]; then rm {1}; fi".format(transientFile)

    return [cmd1, cmd2, cmd3, cmd4]


def removeDuplicates(inputBam, outputBam, cpus=16):
    cmd = "sambamba markdup -t {2} -r {0} {1}".format(inputBam, outputBam, cpus)

    return cmd


def shiftReads(inputBam, genome, outputBam):
    import re

    outputBam = re.sub("\.bam$", "", outputBam)
    cmd1 = "module load samtools"

    cmd2 = "samtools view -h {0} |".format(inputBam)
    cmd2 += " shift_reads.py {0} |".format(genome)
    cmd2 += " samtools view -S -b - |"
    cmd2 += " samtools sort - {0}".format(outputBam)

    return [cmd1, cmd2]


def indexBam(inputBam):
    # module load samtools
    cmd = "samtools index {0}".format(inputBam)

    return cmd


def peakTools(inputBam, output, plot, cpus):
    cmd1 = "module load boost/1.57"
    cmd2 = "module load gcc/4.8.2"
    cmd3 = "module load openmpi/gcc/64/1.8.1"
    cmd4 = "run_spp.R -rf -savp -savp={2} -s=0:5:500 -c={0} -out={1}".format(inputBam, output, plot)

    return [cmd1, cmd2, cmd3, cmd4]


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

    cmd1 = "module load bedtools"

    cmd2 = "bedtools bamtobed -i {0} |".format(inputBam)
    if not tagmented:
        cmd2 += " bedtools slop -i stdin -g {0} -s -l 0 -r 130 |".format(genomeSizes)
        cmd2 += " python fix_bedfile_genome_boundaries.py {0} |".format(genome)
    else:
        cmd2 = "bedtools bamtobed -i {0} |".format(inputBam)
        cmd2 += " get5primePosition.py |"
    cmd2 += " genomeCoverageBed -i stdin -bg -g {0} > {1}.cov".format(genomeSizes, transientFile)

    cmd3 = "bedGraphToBigWig {0}.cov {1} {2}".format(transientFile, genomeSizes, outputBigWig)
    # remove cov file
    cmd4 = "if [[ -s {0}.cov ]]; then rm {0}.cov; fi".format(transientFile)

    cmd5 = "chmod 755 {0}".format(outputBigWig)

    return [cmd1, cmd2, cmd3, cmd4, cmd5]


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


def genomeWideCoverage(inputBam, genomeWindows, output):
    cmd1 = "module load bedtools"
    cmd2 = "bedtools coverage -abam -counts -a {0} -b {1} > {2}".format(inputBam, genomeWindows, output)

    return [cmd1, cmd2]


def calculateFRiP(inputBam, inputBed, output):
    cmd1 = "module load bedtools"

    cmd2 = "cut -f 1,2,3 {0} |".format(inputBed)
    cmd2 += " bedtools coverage -counts -abam {0} -b - |".format(inputBam)
    cmd2 += " awk '{{sum+=$4}} END {{print sum}}' > {0}".format(output)

    return [cmd1, cmd2]


def macs2CallPeaks(treatmentBam, controlBam, outputDir, sampleName, genome, broad=False, plot=True):

    if not broad:
        cmd = "macs2 callpeak -t {0} -c {1} --bw 200 -g {2} -n {3} --outdir {4}""".format(
            treatmentBam, controlBam, genome, sampleName, outputDir
        )
        # --fix-bimodal --extsize 180 \\
    else:
        # Parameter setting for broad factors according to Nature Protocols (2012)
        # Vol.7 No.9 1728-1740 doi:10.1038/nprot.2012.101 Protocol (D) for H3K36me3
        cmd = "macs2 callpeak -t {0} -c {1} --broad --nomodel --extsize 73 --pvalue 1e-3 -g {2} -n {3} --outdir {4}".format(
            treatmentBam, controlBam, genome, sampleName, outputDir
        )

    return cmd


def macs2PlotModel(treatmentBam, controlBam, outputDir, sampleName, genome, broad=False, plot=True):
    import os

    cmd1 = "module load R"

    # run macs r script
    cmd2 = "Rscript {0}/{1}_model.r".format(outputDir, sampleName, os.getcwd())
    # move to sample dir
    cmd3 = "mv {2}/{1}_model.pdf {0}/{1}_model.pdf".format(outputDir, sampleName, os.getcwd())

    return [cmd1, cmd2, cmd3]


def sppCallPeaks(treatmentBam, controlBam, treatmentName, controlName, outputDir, broad, cpus):
    broad = "TRUE" if broad else "FALSE"

    cmd1 = "module load R"

    cmd2 = "Rscript spp_peak_calling.R {0} {1} {2} {3} {4} {5} {6}""".format(
        treatmentBam, controlBam, treatmentName, controlName, broad, cpus, outputDir
    )

    return [cmd1, cmd2]


def bamToBed(inputBam, outputBed):
    cmd1 = "module load bedtools"

    cmd2 = "bedtools bamtobed -i {0} > {1}".format(inputBam, outputBed)

    return [cmd1, cmd2]


def zinbaCallPeaks(treatmentBed, controlBed, cpus, tagmented=False):
    fragmentLength = 80 if tagmented else 180

    cmd1 = "module load R"  # module load R/3.1.0

    cmd2 = "Rscript zinba.R -l {0} -t {1} -c {2}".format(fragmentLength, treatmentBed, controlBed)

    return [cmd1, cmd2]


def homerFindMotifs(peakFile, genome, outputDir, size=150, length="8,10,12,14,16", n_motifs=12):

    cmd = "findMotifsGenome.pl {0} {1} {2}".format(peakFile, genome, outputDir)
    cmd += " -mask -size {3} -len {4} -S {5}".format(size, length, n_motifs)

    return cmd


def AnnotatePeaks(peakFile, genome, motifFile, outputBed):
    cmd = "annotatePeaks.pl {0} {1} -mask -mscore -m {2} |".format(peakFile, genome, motifFile)
    cmd += "tail -n +2 | cut -f 1,5,22 > {3}".format(outputBed)

    return cmd


def centerPeaksOnMotifs(peakFile, genome, windowWidth, motifFile, outputBed):

    cmd = "annotatePeaks.pl {0} {1} -size {2} -center {3} |".format(peakFile, genome, windowWidth, motifFile)
    cmd += " awk -v OFS='\t' '{{print $2, $3, $4, $1, $6, $5}}' |"
    cmd += """ awk -v OFS='\t' -F '\t' '{{ gsub("0", "+", $6) ; gsub("1", "-", $6) ; print }}' |"""
    cmd += " fix_bedfile_genome_boundaries.py {1} | sortBed > {4}".format(genome, outputBed)

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
    cmd += " --window-width {4} --fragment-size {5} --genome {6} --n_clusters {7}".format(
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
    cmd1 = "module load R"

    cmd2 = "Rscript {3}/lib/diffBind_analysis.R {0} {1} {2}".format(
        inputCSV, jobName, plotsDir,
        os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    )

    return [cmd1, cmd2]
