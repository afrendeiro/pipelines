#!/usr/bin/env python

"""
ChIP-seq pipeline 



"""

from argparse import ArgumentParser
import os
import logging
import pandas as pd
import textwrap

__author__ = "Andre Rendeiro"
__copyright__ = "Copyright 2014, Andre Rendeiro"
__credits__ = []
__license__ = "GPL3"
__version__ = "0.1"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"


def slurm_header(jobName, output, queue="shortq", ntasks=1, time="10:00:00", cpusPerTask=16, memPerCpu=2000, nodes=1):
    command = """
    #!/bin/bash
    #SBATCH --partition={queue}
    #SBATCH --ntasks={ntasks}
    #SBATCH --time={time}

    #SBATCH --cpus-per-task={cpusPerTask}
    #SBATCH --mem-per-cpu={memPerTask}
    #SBATCH --nodes={nodes}

    #SBATCH --job-name={jobName}
    #SBATCH --output={output}

    # Start running the job
    hostname
    date

    """.format(queue, ntasks, time, cpusPerTask, memPerCpu, nodes, jobName, output)

    return command


def bam2fastq(inputBam, outputFastq):
    command = """
    # Converting original Bam file to Fastq
    module load bamtools
    
    bamtools convert -in {0} -format fastq > {1}
    
    """.format(inputBam, outputFastq)

    return command


def fastqc(inputBam, outputDir):
    command = """
    # Measuring sample quality with Fastqc
    module load java/jdk/1.7.0_65 
    module load FastQC/0.11.2

    fastqc --noextract --output {0} {1}

    """.format(outputDir, inputBam)

    return command


def trimAdapters(inputFastq, trimmedFastq, adapters):
    command = """
    # Trim adapters from sample
    module load trimmomatic/0.32

    java -jar `which trimmomatic-0.32.jar` SE \
    {0} {1} \
    ILLUMINACLIP:{2}:1:40:15 \
    LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:10 \
    MINLEN:36

    """.format(fastqInput, trimmed, adapters)

    return command

def bowtie2Map(inputFastq, outputBam, genomeIndex, cpus):
    outputBam = re.sub("\.bam" , "", outputBam)
    command = """
    # Map reads with Bowtie2
    module load bowtie/2.2.3
    
    bowtie2 --very-sensitive -p {0} -x {1} {2} | \
    samtools view -S -b - | \
    samtools sort - {3}
    
    """.format(cpus, genomeIndex, inputFastq, outputBam)

    return command


def shiftReads(inputBam, outputBam):
    outputBam = re.sub("\.bam" , "", outputBam)
    command = """
    # Shift read of tagmented sample
    module load bowtie/2.2.3
    
    samtools view -h {0} | \
    python2.7 /home/arendeiro/projects/chipmentation/src/scripts/shift_reads.py | \
    samtools view -S -b - | \
    samtools sort - {1}
    
    """.format(inputBam, outputBam)

    return command
    

def markDuplicates(inputBam, outputBam, metricsFile, tempDir):
    transientFile = re.sub("\.bam" , "", outputBam)
    outputBam = re.sub("\.bam" , "", outputBam)
    command = """
    # Mark duplicates with piccard
    module load bowtie/2.2.3
    
    java -Xmx4g -jar $PICARDDIR/MarkDuplicates.jar
    INPUT={0}
    OUTPUT={1}
    METRICS_FILE={2}
    VALIDATION_STRINGENCY=LENIENT
    TMP_DIR={3}

    # Sort bam file with marked duplicates
    samtools sort {1} {4}

    if [[ -s {4}.bam ]]
        then
        rm {4}.bam
    fi
    
    """.format(inputBam, transientFile, metricsFile, tempDir, outputBam)

    return command


def removeDuplicates():
    pass


def qc():
    $PRESEQ c_curve -B $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
    -o $PROJECTDIR/mapped/$SAMPLE_NAME_qc_c_curve.txt

    $PRESEQ lc_extrap -e 1e8 -s 2e6 -B $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
    -o $PROJECTDIR/mapped/$SAMPLE_NAME_qc_lc_extrap.txt


def mergeReplicates():
    pass


def bamToUCSC():
    pass


def macs2CallPeaks():
    pass


def sppCallPeaks():
    pass


def homerFindMotifs():
    pass


def AnnotatePeaks():
    pass


def centerPeaksOnMotifs():
    pass


def callFootprints():
    pass


def submit_job():
    pass




if __name__ == '__main__':

    # argparser    
    parser = ArgumentParser(
        description = 'ChIP-seq pipeline',
        usage       = 'python chipseq_pipeline.py [OPTIONS] projectName annotation.csv'
        )
    # positional arguments
    parser.add_argument('project_name', help="Project name.", type=str)
    parser.add_argument('csv', help = 'CSV file with sample annotation.', type=str)
    # optional arguments
    parser.add_argument('-r', '--project-root', default = "/fhgfs/groups/lab_bock/shared/projects/",
                        dest = 'project_root', help = 'Specify the directory in which the project will reside.', type=str)
    parser.add_argument('-s', '--stage', default="all", dest='stage',
                        choices=["all", "bam2fastq", "fastqc", "trimming", "mapping",
                                 "shifting", "markduplicates", "removeduplicates", "qc", "mergereplicates"],
                        help='Run only these stages.', type=str)
    parser.add_argument('-g', '--genome', default="hg19", dest='genome',
                        choices=["hg19", "mm10", "zf9"],
                        help="Reference genome.", type=str)
    parser.add_argument('-v', '--verbose', action = 'store_true', dest = 'verbose',
        help = "Verbose behaviour. Print progress to stdout.", default = False)
    
    # parse
    args = parser.parse_args()

    # silent/verbose
    global v
    v = args.verbose

    ### Logging
    # log to current working dir, in the end copy log to projectDir
    logfile_name = os.path.join(os.getcwd(), "log_" + args.project_name + ".log"
    logging.basicConfig(filename=logfile_name,
                        format='%(levelname)s: %(asctime)s %(message)s',
                        datefmt='%Y/%m/%d/ %I:%M:%S %p',
                        level=logging.DEBUG)

    ### Directories and paths
    # check args.project_root exists and user has write access
    args.project_root = os.path.abspath(args.project_root)
    logger.info("Checking if %s directory exists and is writable." % args.project_root)
    if not os.access(args.project_root, os.W_OK):
        log.error("%s does not exist, or user has no write access." % args.project_root)
        sys.exit(1)

    projectDir = os.path.join(args.project_root, args.project_name)
    dataDir = os.path.join(projectDir, "data")

    # make relative project dirs
    dirs = [projectDir,
            dataDir,
            os.path.join(dataDir, "fastq"),
            os.path.join(dataDir, "fastqc"),
            os.path.join(dataDir, "raw"), 
            os.path.join(dataDir, "mapped")
    ]
    for d in dirs:
    if not os.path.exists(d):
        os.makedirs(d)
    
    ### Parse sample information
    # each entry should be in the format:
    # "IP,cellLine,numberCells,technique,experimentName,biologicalReplicate,technicalReplicate"
    args.csv = os.path.abspath(args.csv)
    if not os.path.exists(args.csv):
        log.error("%s does not exist, or user has no write access." % args.csv)
        sys.exit(1)

    samples = pd.read_csv(args.csv)

    # start pipeline
    if args.stage == "all":
        jobCode = slurm_header(jobName, output, queue="shortq", time="10:00:00", cpusPerTask=16, memPerCpu=2000)
        jobCode += bam2fastq(inputBam, outputFastq)
        jobCode += fastqc(inputBam, outputDir)
        jobCode += trimAdapters(fastq, trimmedFastq, "/home/arendeiro/adapters/chipmentation.fa")
        jobCode += bowtie2Map(inputFastq, outputBam, genomeIndex, cpus)
        jobCode += shiftReads(inputBam, outputBam)
        jobCode += markDuplicates(inputBam, outputBam, metricsFile, tempDir)
        jobCode += removeDuplicates()
        jobCode += qc()

    elif args.stage == "bam2fastq":
        jobCode = slurm_header(jobName, output, queue="shortq", time="10:00:00", cpusPerTask=16, memPerCpu=2000)
        jobCode += bam2fastq()
    elif args.stage == "fastqc":
        fastqc()
    elif args.stage == "trimming":
        trimAdapters()
    elif args.stage == "mapping":
        bowtie2Map()
    elif args.stage == "shifting":
        shiftReads()
    elif args.stage == "markduplicates":
        markDuplicates()
    elif args.stage == "removeduplicates":
        removeDuplicates()
    elif args.stage == "qc":
        qc()


    ### if there are samples that are replicates, merge them
    # merge either technical or biological replicates
    # always merge technical before biological - keep all files



    #### DOWNSTREAM

    # Call peaks

    # (Get consensus peaks) or call peaks on merged samples

    # Motif discovery

    # Center peaks on motifs

    # Annotate peaks (various annotations + pwm of motif)

    # Call footprints


    ### SUBMIT
    # Get concatenated string with code from all modules
    jobCode = jobCode

    # Output file name
    jobFile = os.path.join(projectDir, "chipseqPipeline_" + time.strftime("%Y%m%d_%H%M%S") + ".sh")

    with open(jobFile, 'w') as handle:
        handle.write(textwrap.dedent(jobCode))

    # Submit to slurm
    submit_job(jobFile)


    #### 
    # copy log to projectDir
