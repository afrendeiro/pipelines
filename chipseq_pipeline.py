#!/usr/bin/env python

"""
ChIP-seq pipeline 



"""

from argparse import ArgumentParser
import os
import sys
import logging
import pandas as pd
import time
import re
import string
import textwrap
import shutil

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
        logger.error("%s does not exist, or user has no write access." % args.project_root)
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
        os.path.join(resultsDir, "plots")
    ]
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    ### Paths to static files on the cluster
    genomeFolder = "////"
    genomeIndexes = {
        "hg19" : os.path.join(genomeFolder, "hg19"),
        "mm10" : os.path.join(genomeFolder, "mm10"),
        "zf9" : os.path.join(genomeFolder, "zf9")
    }
    adapterFasta = "/home/arendeiro/adapters/chipmentation.fa"
    
    ### Parse sample information
    args.csv = os.path.abspath(args.csv)

    # check it exists
    if not os.path.exists(args.csv):
        logger.error("Sample annotation '%s' does not exist, or user has no read access." % args.csv)
        sys.exit(1)
    
    # read in
    samples = pd.read_csv(args.csv)

    # Perform checks on the variables given
    # (e.g. genome in genomeIndexes)


    # start pipeline
    jobName = string.join([args.project_name, time.strftime("%Y%m%d-%H%M%S")], sep="_")

    # Preprocess each sample individually
    for sample in xrange(len(samples)):
        # if sampleName is not provided, use a concatenation of several variable (excluding longest)
        if "sampleName" in samples.columns and str(samples["sampleName"][sample]) != "nan":
            sampleName = samples["sampleName"][sample]
        else:
            variables = samples.columns.tolist()
            exclude = ["sampleName", "filePath", "genome", "tagmented"]
            [variables.pop(variables.index(exc)) for exc in exclude if exc in variables]
            sampleName = string.join([str(samples[var][sample]) for var in variables], sep="_")
            logger.debug("No sample name provided, using concatenation of variables supplied")

        # check if sample is tagmented or not:
        tagmented = samples["tagmented"][sample] == "yes" or samples["tagmented"][sample] == 1

        # run complete pipeline
        if args.stage == "all":
            jobCode = slurm_header(
                jobName=jobName,
                output=os.path.join(projectDir, "runs", jobName + "_" + sampleName + ".slurm.log"),
                queue="shortq",
                time="10:00:00",
                cpusPerTask=16,
                memPerCpu=2000
            )
        if args.stage in ["all", "bam2fastq"]:
            jobCode += bam2fastq(
                inputBam=samples["filePath"][sample],
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
        if args.stage in ["all", "trimAdapters"]:
            # TODO:
            # Change absolute path to something usable by everyone or to an option.
            jobCode += trimAdapters(
                inputFastq=os.path.join(dataDir, "fastq", sampleName + ".fastq"),
                outputFastq=os.path.join(dataDir, "raw", sampleName + ".trimmed.fastq"),
                adapters=adapterFasta
            )
        if args.stage in ["all", "bowtie2Map"]:
            jobCode += bowtie2Map(
                inputFastq=os.path.join(dataDir, "raw", sampleName + ".trimmed.fastq"),
                outputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.bam"),
                genomeIndex=genomeIndexes[samples["genome"][sample]],
                cpus=args.cpus
            )
        if args.stage in ["all", "shiftReads"]:
            if tagmented:
                # TODO:
                # Get correct relative path.
                # Relative to location of *this* file rather than relative to cwd.
                jobCode += shiftReads(
                    inputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.bam"),
                    outputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.shifted.bam")
                )
        if args.stage in ["all", "markDuplicates"]:
            if tagmented:
                jobCode += markDuplicates(
                    inputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.shifted.bam"),
                    outputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.shifted.dups.bam"),
                    metricsFile=os.path.join(resultsDir, sampleName + ".duplicates.txt")#,
                    #tempDir=
                )
            else:
                jobCode += markDuplicates(
                    inputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.bam"),
                    outputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.dups.bam"),
                    metricsFile=os.path.join(resultsDir, sampleName + ".duplicates.txt")#,
                    #tempDir=
                )
        if args.stage in ["all", "removeDuplicates"]:
            if tagmented:
                jobCode += removeDuplicates(
                    inputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.shifted.dups.bam"),
                    outputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.shifted.nodups.bam")
                )
            else:
                jobCode += removeDuplicates(
                    inputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.dups.bam"),
                    outputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.nodups.bam")
                )
        if args.stage in ["all", "indexBam"]:
            if tagmented:
                jobCode += indexBam(
                    inputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.shifted.dups.bam")
                )
                jobCode += indexBam(
                    inputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.shifted.nodups.bam")
                )
            else:
                jobCode += indexBam(
                    inputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.dups.bam")
                )
                jobCode += indexBam(
                    inputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.nodups.bam")
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
        jobFile = os.path.join(projectDir, "runs", jobName + "_" + sampleName + ".sh")

        with open(jobFile, 'w') as handle:
            handle.write(textwrap.dedent(jobCode))

        # Submit to slurm
        #status = slurm_submit_job(jobFile)

        #if status != 0:
        #    logger.error("Slurm job '%s' not successfull" % jobFile)
        #    sys.exit(1)
        logger.debug("Project '%s'submission finished successfully." % args.project_name)


    ### if there are samples that are replicates, merge them
    # Check if there are replicates

    # merge either technical or biological replicates

    # always merge technical before biological - keep all files


    #### FURTHER TO IMPLEMENT

    # Call peaks

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


def slurm_header(jobName, output, queue="shortq", ntasks=1, time="10:00:00", cpusPerTask=16, memPerCpu=2000, nodes=1):
    command = """    #!/bin/bash
    #SBATCH --partition={0}
    #SBATCH --ntasks={1}
    #SBATCH --time={2}

    #SBATCH --cpus-per-task={3}
    #SBATCH --mem-per-cpu={4}
    #SBATCH --nodes={5}

    #SBATCH --job-name={6}
    #SBATCH --output={7}

    # Start running the job
    hostname
    date

    """.format(queue, ntasks, time, cpusPerTask, memPerCpu, nodes, jobName, output)

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

    java -jar `which trimmomatic-0.32.jar` SE {0} {1}
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
    command = """
    # Shift read of tagmented sample
    module load samtools
    module load python
    
    samtools view -h {0} | \\
    python {1}/lib/shift_reads.py | \\
    samtools view -S -b - | \\
    samtools sort - {2}
    
    """.format(inputBam, os.path.abspath(os.path.curdir), outputBam)

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


def mergeReplicates():
    raise NotImplemented


def bamToUCSC():
    raise NotImplemented


def macs2CallPeaks():
    raise NotImplemented


def sppCallPeaks():
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
    
    # positional arguments
    parser.add_argument('project_name', help="Project name.", type=str)
    parser.add_argument('csv', help='CSV file with sample annotation.', type=str)
    
    # optional arguments
    parser.add_argument('-r', '--project-root', default="/fhgfs/groups/lab_bock/shared/projects/",
                        dest='project_root', type=str,
                        help='Specify the directory in which the project will reside. Default /fhgfs/groups/lab_bock/shared/projects/.')
    parser.add_argument('-s', '--stage', default="all", dest='stage',
                        choices=["all", "bam2fastq", "fastqc", "trimming", "mapping",
                                 "shifting", "markduplicates", "removeduplicates", "qc", "mergereplicates"],
                        help='Run only these stages. Default=all.', type=str)
    parser.add_argument('-c', '--cpus', default=16, dest='cpus',
                        help='Number of CPUs to use. Default=16.', type=int)
    parser.add_argument('-l', '--log-level', default="INFO", dest='log_level',
                        choices=["DEBUG", "INFO", "ERROR"], help='Logging level. Default=INFO.', type=str)
    # parse
    args = parser.parse_args()


    ### Logging
    logger = logging.getLogger(__name__)
    levels = {"DEBUG" : 10, "INFO" : 20, "ERROR" : 40}
    logger.setLevel(levels[args.log_level])

    # create a file handler
    ## (for now in current working dir, in the end copy log to projectDir)
    handler = logging.FileHandler(os.path.join(os.getcwd(), args.project_name + ".log"))
    handler.setLevel(logging.INFO)
    # format logger
    formatter = logging.Formatter(fmt='%(levelname)s: %(asctime)s - %(message)s', datefmt='%Y/%m/%d %H:%M:%S')
    handler.setFormatter(formatter)

    # add the handlers to the logger
    logger.addHandler(handler)


    ### Start main function
    main(args, logger)


    ### Exit
    logger.info("Finished and exiting.")
    sys.exit(0)
