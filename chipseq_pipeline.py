#!/usr/bin/env python

"""
ChIP-seq pipeline 


TODO:
    check if technical replicates (TR), run pipeline for merged TRs
    implement sub command to run either until bams for each biological replicate (BR) or downstream from there
    decide how to treat already existing files (add option, decide default)

    (maybe) implementtt way to check if slurm job is finished and then run 2nd part if there were no errors

"""

from argparse import ArgumentParser
import os
import sys
import logging
import pandas as pd
import numpy as np
import time
import re
import string
import textwrap
import shutil

# TODO: solve pandas chained assignments
pd.options.mode.chained_assignment = None

__author__ = "Andre Rendeiro"
__copyright__ = "Copyright 2014, Andre Rendeiro"
__credits__ = []
__license__ = "GPL3"
__version__ = "0.1"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"

def main(args, logger):
    ### Directories and paths
    # check args.project_root exists and user has write access
    args.project_root = os.path.abspath(args.project_root)
    logger.debug("Checking if %s directory exists and is writable." % args.project_root)
    if not os.access(args.project_root, os.W_OK):
        logger.error("%s does not exist, or user has no write access.\n\
Use option '-r' '--project-root' to set a non-default project root path." % args.project_root)
        sys.exit(1)

        # check args.html_root exists and user has write access
    htmlDir = os.path.abspath(args.html_root)
    logger.debug("Checking if %s directory exists and is writable." % args.project_root)
    if not os.access(htmlDir, os.W_OK):
        logger.error("%s does not exist, or user has no write access.\n\
Use option '--html-root' to set a non-default html root path." % htmlDir)
        sys.exit(1)

    # project directories
    projectDir = os.path.join(args.project_root, args.project_name)
    dataDir = os.path.join(projectDir, "data")
    resultsDir = os.path.join(projectDir, "results")
    
    # make relative project dirs
    dirs = [
        projectDir,
        os.path.join(projectDir, "runs"),
        dataDir,
        os.path.join(dataDir, "fastq"),
        os.path.join(dataDir, "fastqc"),
        os.path.join(dataDir, "raw"), 
        os.path.join(dataDir, "mapped"),
        resultsDir,
        os.path.join(resultsDir, "plots"),
        htmlDir,
        os.path.join(args.html_root, args.project_name),
        os.path.join(args.html_root, args.project_name, "bigWig")
    ]
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)
    htmlDir = os.path.join(args.html_root, args.project_name, "bigWig")
    urlRoot = args.url_root + args.project_name + "/bigWig/"

    ### Paths to static files on the cluster
    genomeFolder = "/fhgfs/prod/ngs_resources/genomes/"
    genomeIndexes = {
        "hg19" : os.path.join(genomeFolder, "hg19/forBowtie2/hg19"),
        "mm10" : os.path.join(genomeFolder, "mm10/forBowtie2/mm10"),
        "zf9" : os.path.join(genomeFolder, "zf9/forBowtie2/zf9")
    }
    genomeSizes = {
        "hg19" : os.path.join(genomeFolder, "hg19/hg19_chromlengths.txt"),
        "mm10" : os.path.join(genomeFolder, "mm10/chrSizes"),
        "zf9" : os.path.join(genomeFolder, "zf9/chrSizes")
    }
    adapterFasta = "/home/arendeiro/adapters/chipmentation.fa"
    
    ### Other static info
    broadFactors = ["H3K27me3", "H3K36me3", "H3K9me3"]


    ### Parse sample information
    args.csv = os.path.abspath(args.csv)

    # check if exists and is a file
    if not os.path.isfile(args.csv):
        logger.error("Sample annotation '%s' does not exist, or user has no read access." % args.csv)
        sys.exit(1)
    
    # read in
    samples = pd.read_csv(args.csv)

    # TODO:
    # Perform checks on the variables given
    # (e.g. genome in genomeIndexes)


    # start pipeline
    projectName = string.join([args.project_name, time.strftime("%Y%m%d-%H%M%S")], sep="_")


    # Check if there are technical replicates
    variables = samples.columns.tolist()
    exclude = ["sampleName", "experimentName", "filePath", "technicalReplicate"]
    [variables.pop(variables.index(exc)) for exc in exclude if exc in variables]

    unique = samples.replace(np.nan, -1).groupby(variables).apply(len).index.values

    biologicalReplicates = dict()
    for sample in xrange(len(unique)):
        replicate = pd.Series(unique[sample],index=variables)
        for row in xrange(len(samples.replace(np.nan, -1)[variables])):
            if (replicate == samples.replace(np.nan, -1)[variables].ix[row]).all():
                if sample not in biologicalReplicates:
                    biologicalReplicates[sample] = [samples.ix[row]]
                else:
                    biologicalReplicates[sample] += [samples.ix[row]]

    # Preprocess biological replicates
    for sample in xrange(len(biologicalReplicates)):
        if len(biologicalReplicates) == 0:
            logger.error("")
            sys.exit(1)

        variables = samples.columns.tolist()
        exclude = ["sampleName", "filePath", "genome", "tagmented"]
        [variables.pop(variables.index(exc)) for exc in exclude if exc in variables]

        # if there's only one technical replicate keep it's name if available
        if len(biologicalReplicates.values()[sample]) == 1:
            s = biologicalReplicates.values()[sample][0]
            # if sampleName is not provided, use a concatenation of several variable (excluding longest)
            if str(s["sampleName"]) != "nan":
                sampleName = samples["sampleName"]
            else:
                sampleName = string.join([str(s[var]) for var in variables], sep="_")
                logger.debug("No sample name provided, using concatenation of variables supplied")
        # if there's more than one technical replicate get a name
        elif len(biologicalReplicates.values()[sample]) > 1:
            s = biologicalReplicates.values()[sample][0]
            s["technicalReplicate"] = 0
            sampleName = string.join([str(s[var]) for var in variables], sep="_")

        jobName = projectName + "_" + sampleName

        # check if sample is tagmented or not:
        tagmented = samples["tagmented"][sample] == "yes" or samples["tagmented"][sample] == 1

        if not tagmented:
            bam = os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.bam")
        else:
            bam = os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.shifted.bam")


        # assemble commands
        # get job header
        jobCode = slurm_header(
            jobName=jobName,
            output=os.path.join(projectDir, "runs", jobName + ".slurm.log"),
            queue=args.queue,
            time=args.time,
            cpusPerTask=args.cpus,
            memPerCpu=args.mem,
            userMail=args.user_mail
        )
        if args.stage in ["all"]:
            # if more than one technical replicate, merge bams
            if len(biologicalReplicates.values()[sample]) > 1:
                jobCode += mergeBams(
                    inputBams=[s["filePath"] for s in biologicalReplicates.values()[sample]],
                    outputBam=os.path.join(dataDir, "raw", sampleName + ".bam")
                )

            # convert bam to fastq
            if args.stage in ["all", "bam2fastq"]:
                if len(biologicalReplicates.values()[sample]) == 1:
                    jobCode += bam2fastq(
                        inputBam=samples["filePath"][sample],
                        outputFastq=os.path.join(dataDir, "fastq", sampleName + ".fastq")
                    )
                elif len(biologicalReplicates.values()[sample]) >= 1:
                    jobCode += bam2fastq(
                        inputBam=os.path.join(dataDir, "raw", sampleName + ".bam"),
                        outputFastq=os.path.join(dataDir, "fastq", sampleName + ".fastq")
                    )
        if args.stage in ["all", "fastqc"]:
            # TODO:
            # Fastqc should be independent from this job but since there's no option in fastqc to specify
            # the sample name, I'll for now run it on the already renamed fastq file produced before,
            # which requires fastqc to run in the same job as the rest :S
            jobCode += fastqc(
                inputFastq=os.path.join(dataDir, "fastq", sampleName + ".fastq"),
                outputDir=os.path.join(dataDir, "fastqc")
            )
        if args.stage in ["all", "trimadapters"]:
            # TODO:
            # Change absolute path to something usable by everyone or to an option.
            jobCode += trimAdapters(
                inputFastq=os.path.join(dataDir, "fastq", sampleName + ".fastq"),
                outputFastq=os.path.join(dataDir, "raw", sampleName + ".trimmed.fastq"),
                adapters=adapterFasta
            )
        if args.stage in ["all", "mapping"]:
            if samples["genome"][sample] not in genomeIndexes:
                logger.error("Sample %s has unsuported genome index: %s" % (sampleName, samples["genome"][sample]))
                sys.exit(1)
            jobCode += bowtie2Map(
                inputFastq=os.path.join(dataDir, "raw", sampleName + ".trimmed.fastq"),
                outputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.bam"),
                genomeIndex=genomeIndexes[samples["genome"][sample]],
                cpus=args.cpus
            )
        if args.stage in ["all", "shiftreads"]:
            if tagmented:
                # TODO:
                # Get correct relative path.
                # Relative to location of *this* file rather than relative to cwd.
                jobCode += shiftReads(
                    inputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.bam"),
                    outputBam=bam
                )
        if args.stage in ["all", "markduplicates"]:
            jobCode += markDuplicates(
                inputBam=bam,
                outputBam=bam + ".dups.bam",
                metricsFile=os.path.join(resultsDir, sampleName + ".duplicates.txt")#,
                #tempDir=
            )
        if args.stage in ["all", "removeduplicates"]:
            jobCode += removeDuplicates(
                inputBam=bam + ".dups.bam",
                outputBam=bam + ".nodups.bam",
            )
        if args.stage in ["all", "indexbam"]:
            jobCode += indexBam(
                inputBam=bam + ".dups.bam"
            )
            jobCode += indexBam(
                inputBam=bam + ".nodups.bam"
            )
        if args.stage in ["all", "maketracks"]:
            jobCode += bamToUCSC(
                inputBam=bam + ".dups.bam",
                outputBigWig=os.path.join(htmlDir, sampleName + ".bigWig"),
                genomeSizes=genomeSizes[samples["genome"][sample]],
                tagmented=False
            )
            jobCode += addTrackToHub(
                sampleName=sampleName,
                trackURL=urlRoot + sampleName + ".bigWig",
                trackHub=os.path.join(htmlDir, "trackHub.txt")
            )
            jobCode += bamToUCSC(
                inputBam=bam + ".dups.bam",
                outputBigWig=os.path.join(htmlDir, sampleName + ".5prime.bigWig"),
                genomeSizes=genomeSizes[samples["genome"][sample]],
                tagmented=True
            )
            jobCode += addTrackToHub(
                sampleName=sampleName,
                trackURL=urlRoot + sampleName + ".5prime.bigWig",
                trackHub=os.path.join(htmlDir, "trackHub.txt")
            )

        # if args.stage in ["all", "qc"]:
        #     if tagmented:
        #         jobCode += qc()
        #     else:
        #         jobCode += qc()

        ### Submit job to slurm
        # Get concatenated string with code from all modules
        jobCode += slurm_footer()

        # Output file name
        jobFile = os.path.join(projectDir, "runs", jobName + ".sh")

        with open(jobFile, 'w') as handle:
            handle.write(textwrap.dedent(jobCode))

        # Submit to slurm
        #status = slurm_submit_job(jobFile)

        #if status != 0:
        #    logger.error("Slurm job '%s' not successfull" % jobFile)
        #    sys.exit(1)
        logger.debug("Project '%s'submission finished successfully." % args.project_name)


    #### FURTHER TO IMPLEMENT

    # Call peaks
    # if args.stage in ["all", "callPeaks"]:
    #     jobCode += macs2CallPeaks(
    #         inputBam=os.path.join(
    #         treatmentBam=bam + "dups.bam",
    #         controlBam=,
    #         outputDir=,
    #         sampleName=,
    #         broad=False
    #     )

    # (Get consensus peaks) or call peaks on merged samples

    # Motif discovery

    # Center peaks on motifs

    # Annotate peaks (various annotations + pwm of motif)

    # Call footprints
    

    #### Copy log to projectDir
    shutil.copy2(
        os.path.join(os.getcwd(), args.project_name + ".log"),
        os.path.join(projectDir, "runs", args.project_name + ".log")
    )
    logger.debug("Copied log file to project directory '%s'" % os.path.join(projectDir, "runs"))


def slurm_header(jobName, output, queue="shortq", ntasks=1, time="10:00:00", cpusPerTask=16, memPerCpu=2000, nodes=1, userMail=""):
    command = """    #!/bin/bash
    #SBATCH --partition={0}
    #SBATCH --ntasks={1}
    #SBATCH --time={2}

    #SBATCH --cpus-per-task={3}
    #SBATCH --mem-per-cpu={4}
    #SBATCH --nodes={5}

    #SBATCH --job-name={6}
    #SBATCH --output={7}

    #SBATCH --mail-type=fail,end
    #SBATCH --mail-user={8}

    # Start running the job
    hostname
    date

    """.format(queue, ntasks, time, cpusPerTask, memPerCpu, nodes, jobName, output, userMail)

    return command


def slurm_footer():
    command = """
    # Job end
    date

    """

    return command


def slurm_submit_job(jobFile):
    command = "sbatch %s" % jobFile

    return os.system(command)


def mergeBams(inputBams, outputBam):
    command = """
    # Merging bam files from replicates
    module load samtools

    samtools merge {0} {1}

    """.format(outputBam, (" ".join(["%s"] * len(inputBams))) % tuple(inputBams))

    return command


def bam2fastq(inputBam, outputFastq):
    command = """
    # Converting original Bam file to Fastq
    module load bamtools
    
    bamtools convert -in {0} -format fastq -out {1}
    
    """.format(inputBam, outputFastq)

    return command


def fastqc(inputFastq, outputDir):
    command = """
    # Measuring sample quality with Fastqc
    module load java/jdk/1.7.0_65 
    module load FastQC/0.11.2

    fastqc --noextract --output {0} {1}

    """.format(outputDir, inputFastq)

    return command


def trimAdapters(inputFastq, outputFastq, adapters):
    command = """
    # Trim adapters from sample
    module load trimmomatic/0.32

    java -jar `which trimmomatic-0.32.jar` SE {0} {1} \\
    ILLUMINACLIP:{2}:1:40:15 \\
    LEADING:3 TRAILING:3 \\
    SLIDINGWINDOW:4:10 \\
    MINLEN:36

    """.format(inputFastq, outputFastq, adapters)

    return command

def bowtie2Map(inputFastq, outputBam, genomeIndex, cpus):
    outputBam = re.sub("\.bam^" , "", outputBam)
    command = """
    # Map reads with Bowtie2
    module load bowtie/2.2.3
    
    bowtie2 --very-sensitive -p {0} -x {1} {2} | \\
    samtools view -S -b - | \\
    samtools sort - {3}
    
    """.format(cpus, genomeIndex, inputFastq, outputBam)

    return command


def shiftReads(inputBam, outputBam):
    outputBam = re.sub("\.bam" , "", outputBam)
    ### TODO:
    # Implement read shifting with HTSeq or Cython
    command = """
    # Shift read of tagmented sample
    module load samtools
    module load python
    
    samtools view -h {0} | \\
    python {1}/lib/shift_reads.py | \\
    samtools view -S -b - | \\
    samtools sort - {2}
    
    """.format(inputBam, os.path.abspath(os.path.dirname(os.path.realpath(__file__))), outputBam)

    return command
    

def markDuplicates(inputBam, outputBam, metricsFile, tempDir="."):
    transientFile = re.sub("\.bam" , "", outputBam)
    outputBam = re.sub("\.bam" , "", outputBam)
    command = """
    # Mark duplicates with piccard
    module load bowtie/2.2.3
    
    java -Xmx4g -jar $PICARDDIR/MarkDuplicates.jar \\
    INPUT={0} \\
    OUTPUT={1} \\
    METRICS_FILE={2} \\
    VALIDATION_STRINGENCY=LENIENT \\
    TMP_DIR={3}

    # Sort bam file with marked duplicates
    samtools sort {1} {4}

    if [[ -s {4}.bam ]]
        then
        rm {4}.bam
    fi
    
    """.format(inputBam, transientFile, metricsFile, tempDir, outputBam)

    return command


def removeDuplicates(inputBam, outputBam):
    command = """
    # Remove duplicates with sambamba
    sambamba markdup -t 16 -r {0} {1}

    """.format(inputBam, outputBam)

    return command


def indexBam(inputBam):
    command = """
    # Indexing bamfile with samtools
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
    raise NotImplemented


def bamToUCSC(inputBam, outputBigWig, genomeSizes, tagmented=False):
    transientFile = os.path.abspath(re.sub("\.bigWig" , "", outputBigWig))
    if not tagmented:
        command = """
    # make bigWig tracks from bam file
    module load 

    bamToBed -i {0} | \\
    genomeCoverageBed -i stdin -bg -g {1} > {2}.cov

    bedGraphToBigWig {2}.cov {1} {3}

    # remove cov file
    if [[ -s {2}.cov ]]
        then
        rm {2}.cov
    fi

    """.format(inputBam, genomeSizes, transientFile, outputBigWig)
    else:
        command = """
    # make bigWig tracks from bam file
    module load 

    bamToBed -i {0} | \\
    python {4}/lib/get5primePosition.py | \\
    genomeCoverageBed -i stdin -bg -g {1} > {2}.cov

    bedGraphToBigWig {2}.cov {1} {3}

    # remove cov file
    if [[ -s {2}.cov ]]
        then
        rm {2}.cov
    fi

    """.format(inputBam, genomeSizes, transientFile, outputBigWig, os.path.abspath(os.path.dirname(os.path.realpath(__file__))))

    return command


def addTrackToHub(sampleName, trackURL, trackHub):
    command = """
    echo 'track type=bigWig name="{0}" description="{0}" visibility=3 bigDataUrl={1}' >> {2}
    """.format(sampleName, trackURL, trackHub)

    return command


def macs2CallPeaks(treatmentBam, controlBam, outputDir, sampleName, broad=False):
    if not broad:
        command = """
    # call peaks with MACS2

    macs2 callpeak -t {0} -c {1} \\
    --bw 200 \\
    -g hs -n {2} --outdir {3}

    """.format(treatmentBam, controlBam, sampleName, outputDir)
    else:
        command = """
    # call peaks with MACS2
        
    macs2 callpeak -B -t {0} -c {1} \\
    --broad --nomodel --extsize 200 --pvalue 1e-3 \\
    -g hs -n {2} --outdir {3}

    """.format(treatmentBam, controlBam, sampleName, outputDir)

    return command


def sppCallPeaks(treatmentBam, controlBam, outputDir, broad=False):
    raise NotImplemented


def homerFindMotifs():
    raise NotImplemented


def AnnotatePeaks():
    raise NotImplemented


def centerPeaksOnMotifs():
    raise NotImplemented


def callFootprints():
    raise NotImplemented



if __name__ == '__main__':
    ### Parse command-line arguments
    parser = ArgumentParser(
        description = 'ChIP-seq pipeline',
        usage       = 'python chipseq_pipeline.py [OPTIONS] projectName annotation.csv'
    )

    # subcommands
    # subparsers = parser.add_subparsers(title='Sub commands',
    #                                description='Valid subcommands',
    #                                help='individual: Process biological replicate samples individualy. comparison: Merge biological replicates and perform comparisons.')
    # subparsers.add_parser('individual')
    # subparsers.add_parser('comparison')
    # parser.parse_args(['-h'])
    
    # positional arguments
    parser.add_argument('project_name', help="Project name.", type=str)
    parser.add_argument('csv', help='CSV file with sample annotation.', type=str)
    
    # optional arguments
    parser.add_argument('-r', '--project-root', default="/fhgfs/groups/lab_bock/shared/projects/",
                        dest='project_root', type=str,
                        help='Directory in which the project will reside. Default=/fhgfs/groups/lab_bock/shared/projects/.')
    parser.add_argument('-s', '--stage', default="all", dest='stage',
                        choices=["all", "bam2fastq", "fastqc", "trimming", "mapping",
                                 "shiftreads", "markduplicates", "removeduplicates", "indexbam", "qc", "maketracks", "mergereplicates"],
                        help='Run only these stages. Default=all.', type=str)
    parser.add_argument('-c', '--cpus', default=16, dest='cpus',
                        help='Number of CPUs to use. Default=16.', type=int)
    parser.add_argument('-m', '--mem-per-cpu', default=2000, dest='mem',
                        help='Memory per CPU to use. Default=2000.', type=int)
    parser.add_argument('-q', '--queue', default="shortq", dest='queue',
                        choices=["shortq", "mediumq", "longq"],
                        help='Queue to submit jobs to. Default=shortq', type=str)
    parser.add_argument('-t', '--time', default="10:00:00", dest='time',
                        help='Maximum time for jobs to run. Default=10:00:00', type=str)
    parser.add_argument('--user-mail', default="", dest='user_mail',
                        help='User mail address. Default=<submitting user>.', type=str)
    parser.add_argument('--html-root', default="/fhgfs/groups/lab_bock/public_html/",
                        dest='html_root', type=str,
                        help='public_html directory in which bigwig files for the project will reside. Default=/fhgfs/groups/lab_bock/public_html/.')
    parser.add_argument('--url-root', default="http://www.biomedical-sequencing.at/bocklab/",
                        dest='url_root', type=str,
                        help='Url mapping to public_html directory where bigwig files for the project will be accessed. Default=http://www.biomedical-sequencing.at/bocklab/.')
    parser.add_argument('-l', '--log-level', default="INFO", dest='log_level',
                        choices=["DEBUG", "INFO", "ERROR"], help='Logging level. Default=INFO.', type=str)
    # parse
    args = parser.parse_args()


    ### Logging
    logger = logging.getLogger(__name__)
    levels = {"DEBUG" : 10, "INFO" : 20, "ERROR" : 40}
    logger.setLevel(levels[args.log_level])

    # create a file handler
    # (for now in current working dir, in the end copy log to projectDir)
    handler = logging.FileHandler(os.path.join(os.getcwd(), args.project_name + ".log"))
    handler.setLevel(logging.INFO)
    # format logger
    formatter = logging.Formatter(fmt='%(levelname)s: %(asctime)s - %(message)s', datefmt='%Y/%m/%d %H:%M:%S')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # create a stout handler
    stdout = logging.StreamHandler(sys.stdout)
    stdout.setLevel(logging.ERROR)
    formatter = logging.Formatter(fmt='%(levelname)s: %(asctime)s - %(message)s', datefmt='%Y/%m/%d %H:%M:%S')
    stdout.setFormatter(formatter)
    logger.addHandler(stdout)

    ### Start main function
    main(args, logger)


    ### Exit
    logger.info("Finished and exiting.")
    sys.exit(0)
