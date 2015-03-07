#!/usr/bin/env python

"""
ChIP-seq pipeline

    ## BUGS
    Take care of additional requirements: UCSCutils, ...

    ## FURTHER TO IMPLEMENT
    copy original sample annotation file to projectDir
    merge biological replicates, process again
    Rscript macs2_model.r

"""

from argparse import ArgumentParser
import os
import sys
import logging
import pandas as pd
import numpy as np
import time
import random
import re
import string
import textwrap
import shutil
import subprocess

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


def main():
    # Parse command-line arguments
    parser = ArgumentParser(description='ChIP-seq pipeline.')

    # Global options
    # positional arguments
    # optional arguments
    parser.add_argument('-r', '--project-root', default="/fhgfs/groups/lab_bock/shared/projects/",
                        dest='project_root', type=str,
                        help='Directory in which the project will reside. Default=/fhgfs/groups/lab_bock/shared/projects/.')
    parser.add_argument('--html-root', default="/fhgfs/groups/lab_bock/public_html/arendeiro/",
                        dest='html_root', type=str,
                        help='public_html directory in which bigwig files for the project will reside. Default=/fhgfs/groups/lab_bock/public_html/.')
    parser.add_argument('--url-root', default="http://www.biomedical-sequencing.at/bocklab/arendeiro/",
                        dest='url_root', type=str,
                        help='Url mapping to public_html directory where bigwig files for the project will be accessed. Default=http://www.biomedical-sequencing.at/bocklab.')
    parser.add_argument('--keep-tmp-files', dest='keep_tmp', action='store_true',
                        help='Keep intermediary files. If not it will only preserve final files. Default=False.')
    parser.add_argument('-c', '--cpus', default=16, dest='cpus',
                        help='Number of CPUs to use. Default=16.', type=int)
    parser.add_argument('-m', '--mem-per-cpu', default=2000, dest='mem',
                        help='Memory per CPU to use. Default=2000.', type=int)
    parser.add_argument('-q', '--queue', default="shortq", dest='queue',
                        choices=["develop", "shortq", "mediumq", "longq"],
                        help='Queue to submit jobs to. Default=shortq', type=str)
    parser.add_argument('-t', '--time', default="10:00:00", dest='time',
                        help='Maximum time for jobs to run. Default=10:00:00', type=str)
    parser.add_argument('--user-mail', default="", dest='user_mail',
                        help='User mail address. Default=<submitting user>.', type=str)
    parser.add_argument('-l', '--log-level', default="INFO", dest='log_level',
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help='Logging level. Default=INFO.', type=str)
    parser.add_argument('--dry-run', dest='dry_run', action='store_true',
                        help='Dry run. Assemble commands, but do not submit jobs to slurm. Default=False.')

    # Sub commands
    subparser = parser.add_subparsers(title='sub-command', dest="command")
    # preprocess
    preprocess_subparser = subparser.add_parser("preprocess")
    preprocess_subparser.add_argument(dest='project_name', help="Project name.", type=str)
    preprocess_subparser.add_argument(dest='csv', help='CSV file with sample annotation.', type=str)
    preprocess_subparser.add_argument('-s', '--stage', default="all", dest='stage',
                                      choices=["all", "bam2fastq", "fastqc", "trimming", "mapping",
                                               "shiftreads", "markduplicates", "removeduplicates",
                                               "indexbam", "qc", "maketracks"],
                                      help='Run only these stages. Default=all.', type=str)
    # analyse
    comparison_subparser = subparser.add_parser("analyse")
    comparison_subparser.add_argument(dest='project_name', help="Project name.", type=str)
    comparison_subparser.add_argument(dest='csv', help='CSV file with sample annotation.', type=str)
    comparison_subparser.add_argument('-s', '--stage', default="all", dest='stage',
                                      choices=["all", "callpeaks", "findmotifs", "centerpeaks",
                                               "annotatepeaks", "peakanalysis", "tssanalysis", "footprints"],
                                      help='Run only these stages. Default=all.', type=str)
    comparison_subparser.add_argument('--peak-caller', default="macs2", choices=["macs2", "spp"],
                                      dest='peak_caller', help='Peak caller to use. Default=macs2.', type=str)
    comparison_subparser.add_argument('--peak-window-width', default=2000,
                                      dest='peak_window_width', help='Width of window around peak motifs. Default=2000.', type=int)
    comparison_subparser.add_argument('--duplicates', dest='duplicates', action='store_true',
                                      help='Allow duplicates in coorelation analysis. Default=False.')
    comparison_subparser.add_argument('--genome-window-width', default=1000,
                                      dest='genome_window_width', help='Width of window to make genome-wide correlations. Default=1000.', type=int)

    # compare
    comparison_subparser = subparser.add_parser("compare")
    comparison_subparser.add_argument(dest='project_name', help="Project name.", type=str)
    comparison_subparser.add_argument(dest='csv', help='CSV file with sample annotation.', type=str)
    comparison_subparser.add_argument('-s', '--stage', default="all", dest='stage',
                                      choices=["all", "diffbind", "correlations"],
                                      help='Run only these stages. Default=all.', type=str)

    # Parse
    args = parser.parse_args()

    # Logging
    logger = logging.getLogger(__name__)
    levels = {"DEBUG": logging.DEBUG, "INFO": logging.INFO, 'WARNING': logging.WARNING, "ERROR": logging.ERROR, "CRITICAL": logging.CRITICAL}
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

    # Start main function
    if args.command == "preprocess":
        logger, fr, to = preprocess(args, logger)
    elif args.command == "stats":
        logger, fr, to = readStats(args, logger)
    elif args.command == "analyse":
        logger, fr, to = analyse(args, logger)

    elif args.command == "compare":
        logger, fr, to = compare(args, logger)

    # Exit
    logger.info("Finished and exiting.")

    # Copy log to projectDir
    shutil.copy2(fr, to)

    sys.exit(1)


def checkProjectDirs(args, logger):
    # Directories and paths
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
        os.path.join(dataDir, "coverage"),
        os.path.join(dataDir, "peaks"),
        os.path.join(dataDir, "motifs"),
        resultsDir,
        os.path.join(resultsDir, "plots"),
        htmlDir,
        os.path.join(args.html_root, args.project_name),
        os.path.join(args.html_root, args.project_name, "bigWig")
    ]
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    # chmod of paths to public_html folder
    html = [
        htmlDir,
        os.path.join(args.html_root, args.project_name),
        os.path.join(args.html_root, args.project_name, "bigWig")
    ]
    for d in html:
        try:
            os.chmod(d, 0755)
        except OSError:
            logger.error("cannot change folder's mode: %s" % d)
            continue

    htmlDir = os.path.join(args.html_root, args.project_name, "bigWig")
    urlRoot = args.url_root + args.project_name + "/bigWig/"

    return (htmlDir, projectDir, dataDir, resultsDir, urlRoot)


def checkTechnicalReplicates(samples):
    """
    returns dictionary with entries and list of pd.Series with sample info
    samples - a pandas.DataFrame with sample info.

    """
    # Check if there are technical replicates
    variables = samples.columns.tolist()
    exclude = ["sampleNumber", "sampleName", "technicalReplicate", "filePath", "controlSampleName"]
    [variables.pop(variables.index(exc)) for exc in exclude if exc in variables]

    unique = samples.replace(np.nan, -1).groupby(variables).apply(len).index.values

    biologicalReplicates = dict()
    for sample in xrange(len(unique)):
        replicate = pd.Series(unique[sample], index=variables)
        for row in xrange(len(samples.replace(np.nan, -1)[variables])):
            if (replicate == samples.replace(np.nan, -1)[variables].ix[row]).all():
                if sample not in biologicalReplicates:
                    biologicalReplicates[sample] = [samples.ix[row]]
                else:
                    biologicalReplicates[sample] += [samples.ix[row]]

    return biologicalReplicates


def preprocess(args, logger):
    logger.info("Starting sample preprocessing.")

    logger.debug("Checking project directories exist and creating if not.")
    htmlDir, projectDir, dataDir, resultsDir, urlRoot = checkProjectDirs(args, logger)

    # Paths to static files on the cluster
    genomeFolder = "/fhgfs/prod/ngs_resources/genomes/"
    genomeIndexes = {
        "hg19": os.path.join(genomeFolder, "hg19/forBowtie2/withoutRandom/hg19"),
        "mm10": os.path.join(genomeFolder, "mm10/forBowtie2/mm10/index_woRandom/mm10"),
        "dr7": os.path.join(genomeFolder, "dr7/forBowtie2/dr7")
    }
    genomeSizes = {
        "hg19": "/fhgfs/groups/lab_bock/arendeiro/share/hg19.chrom.sizes",
        "mm10": "/fhgfs/groups/lab_bock/arendeiro/share/mm10.chrom.sizes",
        "dr7": "/fhgfs/groups/lab_bock/arendeiro/share/danRer7.chrom.sizes"
    }
    adapterFasta = "/fhgfs/groups/lab_bock/shared/cm.fa"

    # Other static info
    tagment = [
        "DNASE", "DNASESEQ", "DNASE-SEQ", "DHS", "DHS-SEQ", "DHSSEQ",
        "ATAC", "ATAC-SEQ", "ATACSEQ",
        "CM"
    ]

    # from the ggplot2 color blind pallete
    # #999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
    # shortly: green,red=active, orange=transcription, blue=repression, yellow=potential
    colours = {
        "IGG": "153,153,153", "INPUT": "153,153,153",  # grey
        "H3K36ME1": "230,159,0", "H3K36ME2": "230,159,0", "H3K36ME3": "230,159,0",  # orange
        "H3K4ME3": "0,158,115",  # bluish green
        "H3K4ME1": "120,114,33", "H3K14ac": "120,114,33",  # yellow
        "H3K27ME1": "0,114,178", "H3K27ME2": "0,114,178", "H3K27ME3": "0,114,178",  # blue
        "H3K9ME1": "86,180,233", "H3K9ME2": "86,180,233", "H3K9ME3": "86,180,233",  # sky blue
        "H3AC": "213,94,0", "H3K9AC": "213,94,0", "H3K27AC": "213,94,0", "H3K56AC": "213,94,0", "H3K56AC": "213,94,0",  # vermillion
        "H3K79ME1": "204,121,167", "H3K79ME2": "204,121,167", "H3K79ME3": "204,121,167"  # reddish purple
    }

    colour_gradient = [  # 10 colour gradient from red to blue
        "155,3,5", "140,2,18", "125,2,31", "110,2,44", "96,2,57",
        "81,2,70", "66,2,83", "52,2,96", "37,2,109", "22,2,122"
    ]

    # Parse sample information
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

    # Get biological replicates from technical
    logger.debug("Checking which technical replicates form biological replicates.")
    biologicalReplicates = checkTechnicalReplicates(samples)  # sorted differently than input annotation

    # Preprocess biological replicates
    samplesMerged = list()
    for sample in xrange(len(biologicalReplicates)):
        if len(biologicalReplicates) == 0:
            logger.error("No samples in sheet.")
            sys.exit(1)

        # Get sample name
        variables = samples.columns.tolist()
        exclude = ["sampleNumber", "sampleName", "filePath", "genome", "controlSampleName"]
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
            s = biologicalReplicates[sample][0]
            s["technicalReplicate"] = 0
            sampleName = string.join([str(s[var]) for var in variables], sep="_")

        # add sample name to series
        s = biologicalReplicates[sample][0]
        s['sampleName'] = sampleName

        # get jobname
        jobName = projectName + "_" + sampleName

        # check if sample is paired end
        PE = checkSamplePE(biologicalReplicates.values()[0]["filePath"])

        # check if sample is tagmented or not:
        tagmented = True if biologicalReplicates[sample][0]["technique"] in tagment else False

        # get intermediate names for files
        if not tagmented:
            bam = os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2")
        else:
            bam = os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.shifted")

        if len(biologicalReplicates.values()[sample]) == 1:
            unmappedBam = biologicalReplicates.values()[sample][0]["filePath"]
        else:
            unmappedBam = os.path.join(dataDir, "raw", sampleName + ".bam")

        # get colour for tracks
        if samples["ip"][sample].upper() in colours.keys():
            colour = colours[samples["ip"][sample].upper()]
        else:
            colour = random.sample(colour_gradient, 1)[0]  # pick one randomly

        # append to list of Biological Replicates
        toAppend = s.copy()
        toAppend['filePath'] = bam + ".dups.bam"  # bam file of merged sample if several technical replicates
        samplesMerged.append(toAppend)

        # keep track of temporary files
        tempFiles = list()

        # assemble commands
        # get job header
        jobCode = slurmHeader(
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
                    inputBams=[ss["filePath"] for ss in biologicalReplicates.values()[sample]],
                    outputBam=unmappedBam
                )
        if args.stage in ["all", "fastqc"]:
            # TODO:
            # Fastqc should be independent from this job but since there's no option in fastqc to specify
            # the sample name, I'll for now run it on the already renamed fastq file produced before,
            # which requires fastqc to run in the same job as the rest :S
            jobCode += fastqc(
                inputBam=unmappedBam,
                outputDir=os.path.join(dataDir, "fastqc")
            )
        # convert bam to fastq
        if args.stage in ["all", "bam2fastq"]:
            if not PE:
                jobCode += bam2fastq(
                    inputBam=unmappedBam,
                    outputFastq=os.path.join(dataDir, "fastq", sampleName + ".fastq")
                )
                tempFiles.append(os.path.join(dataDir, "fastq", sampleName + ".fastq"))
            else:
                jobCode += bam2fastq(
                    inputBam=unmappedBam,
                    outputFastq=os.path.join(dataDir, "fastq", sampleName + "_1.fastq"),
                    outputFastq2=os.path.join(dataDir, "fastq", sampleName + "_2.fastq"),
                    unpairedFastq=os.path.join(dataDir, "fastq", sampleName + "_unpaired.fastq")
                )
                tempFiles.append(os.path.join(dataDir, "fastq", sampleName + "_1.fastq"))
                tempFiles.append(os.path.join(dataDir, "fastq", sampleName + "_2.fastq"))
                tempFiles.append(os.path.join(dataDir, "fastq", sampleName + "_unpaired.fastq"))

        if args.stage in ["all", "trimadapters"]:
            # TODO:
            # Change absolute path to something usable by everyone or to an option.
            if not PE:
                jobCode += trimAdaptersSE(
                    inputFastq=os.path.join(dataDir, "fastq", sampleName + ".fastq"),
                    outputFastq=os.path.join(dataDir, "fastq", sampleName + ".trimmed.fastq"),
                    cpus=args.cpus,
                    adapters=adapterFasta,
                    log=os.path.join(resultsDir, sampleName + ".trimlog")
                )
                tempFiles.append(os.path.join(dataDir, "fastq", sampleName + ".trimmed.fastq"))
            else:
                jobCode += trimAdaptersPE(
                    inputFastq1=os.path.join(dataDir, "fastq", sampleName + "_1.fastq"),
                    inputFastq2=os.path.join(dataDir, "fastq", sampleName + "_2.fastq"),
                    outputFastq1=os.path.join(dataDir, "fastq", sampleName + "_1.trimmed.fastq"),
                    outputFastq1unpaired=os.path.join(dataDir, "fastq", sampleName + "_1_unpaired.trimmed.fastq"),
                    outputFastq2=os.path.join(dataDir, "fastq", sampleName + "_2.trimmed.fastq"),
                    outputFastq2unpaired=os.path.join(dataDir, "fastq", sampleName + "_2_unpaired.trimmed.fastq"),
                    cpus=args.cpus,
                    adapters=adapterFasta,
                    log=os.path.join(resultsDir, sampleName + ".trimlog")
                )
                tempFiles.append(os.path.join(dataDir, "fastq", sampleName + "_1.trimmed.fastq"))
                tempFiles.append(os.path.join(dataDir, "fastq", sampleName + "_2.trimmed.fastq"))
                tempFiles.append(os.path.join(dataDir, "fastq", sampleName + "_1_unpaired.trimmed.fastq"))
                tempFiles.append(os.path.join(dataDir, "fastq", sampleName + "_2_unpaired.trimmed.fastq"))

        if args.stage in ["all", "mapping"]:
            if samples["genome"][sample] not in genomeIndexes:
                logger.error("Sample %s has unsuported genome index: %s" % (sampleName, samples["genome"][sample]))
                sys.exit(1)
            if not PE:
                jobCode += bowtie2Map(
                    inputFastq=os.path.join(dataDir, "fastq", sampleName + ".trimmed.fastq"),
                    outputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.bam"),
                    genomeIndex=genomeIndexes[samples["genome"][sample]],
                    cpus=args.cpus
                )
            else:
                jobCode += bowtie2Map(
                    inputFastq=os.path.join(dataDir, "fastq", sampleName + "_1.trimmed.fastq"),
                    outputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.bam"),
                    genomeIndex=genomeIndexes[samples["genome"][sample]],
                    cpus=args.cpus,
                    inputFastq2=os.path.join(dataDir, "fastq", sampleName + "_2.trimmed.fastq")
                )
            tempFiles.append(os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.bam"))
        if args.stage in ["all", "shiftreads"]:
            if tagmented:
                jobCode += shiftReads(
                    inputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.bam"),
                    genome=samples["genome"][sample],
                    outputBam=bam + ".bam"
                )
        if args.stage in ["all", "markduplicates"]:
            jobCode += markDuplicates(
                inputBam=bam + ".bam",
                outputBam=bam + ".dups.bam",
                metricsFile=os.path.join(resultsDir, sampleName + ".duplicates.txt")  #
                # tempDir=
            )
        if args.stage in ["all", "removeduplicates"]:
            jobCode += removeDuplicates(
                inputBam=bam + ".dups.bam",
                outputBam=bam + ".nodups.bam",
                cpus=args.cpus
            )
        if args.stage in ["all", "indexbam"]:
            jobCode += indexBam(
                inputBam=bam + ".dups.bam"
            )
            jobCode += indexBam(
                inputBam=bam + ".nodups.bam"
            )
        if args.stage in ["all", "maketracks"]:
            # right now tracks are only made for bams with duplicates
            jobCode += bamToUCSC(
                inputBam=bam + ".dups.bam",
                outputBigWig=os.path.join(htmlDir, sampleName + ".bigWig"),
                genomeSizes=genomeSizes[samples["genome"][sample]],
                genome=samples["genome"][sample],
                tagmented=False
            )
            jobCode += addTrackToHub(
                sampleName=sampleName,
                trackURL=urlRoot + sampleName + ".bigWig",
                trackHub=os.path.join(htmlDir, "trackHub_{0}.txt".format(samples["genome"][sample])),
                colour=colour
            )
            if tagmented:
                jobCode += bamToUCSC(
                    inputBam=bam + ".dups.bam",
                    outputBigWig=os.path.join(htmlDir, sampleName + ".5prime.bigWig"),
                    genomeSizes=genomeSizes[samples["genome"][sample]],
                    genome=samples["genome"][sample],
                    tagmented=True
                )
                jobCode += addTrackToHub(
                    sampleName=sampleName,
                    trackURL=urlRoot + sampleName + ".5prime.bigWig",
                    trackHub=os.path.join(htmlDir, "trackHub_{0}.txt".format(samples["genome"][sample])),
                    fivePrime="5prime",
                    colour=colour
                )
            linkToTrackHub(
                trackHubURL="{0}/{1}/bigWig/trackHub_{2}.txt".format(args.url_root, args.project_name, samples["genome"][sample]),
                fileName=os.path.join(projectDir, "ucsc_tracks_{0}.html".format(samples["genome"][sample])),
                genome=samples["genome"][sample]
            )

        # if args.stage in ["all", "qc"]:
        #     if tagmented:
        #         jobCode += qc()
        #     else:
        #         jobCode += qc()

        # Remove intermediary files
        if args.stage == "all" and not args.keep_tmp:
            logger.debug("Removing intermediary files")
            for fileName in tempFiles:
                jobCode += removeFile(fileName)

        # Submit job to slurm
        # Get concatenated string with code from all modules
        jobCode += slurmFooter()

        # Output file name
        jobFile = os.path.join(projectDir, "runs", jobName + ".sh")

        with open(jobFile, 'w') as handle:
            handle.write(textwrap.dedent(jobCode))

        # Submit to slurm
        if not args.dry_run:
            logger.info("Submitting jobs to slurm")
            status = slurmSubmitJob(jobFile)

            if status != 0:
                logger.error("Slurm job '%s' not successfull" % jobFile)
                sys.exit(1)
            logger.debug("Project '%s'submission finished successfully." % args.project_name)

    # write annotation sheet with biological replicates
    df = pd.DataFrame(samplesMerged)
    df["controlSampleName"] = None
    df.to_csv(os.path.join(projectDir, args.project_name + ".biol_replicates.annotation_sheet.csv"), index=False)

    logger.debug("Finished preprocessing")

    return (logger,
            os.path.join(os.getcwd(), args.project_name + ".log"),
            os.path.join(projectDir, "runs", args.project_name + ".log"))


def readStats(args, logger):
    logger.info("Starting sample read stats.")

    logger.debug("Checking project directories exist and creating if not.")
    htmlDir, projectDir, dataDir, resultsDir, urlRoot = checkProjectDirs(args, logger)

    # Parse sample information
    args.csv = os.path.abspath(args.csv)

    # check if exists and is a file
    if not os.path.isfile(args.csv):
        logger.error("Sample annotation '%s' does not exist, or user has no read access." % args.csv)
        sys.exit(1)

    # read in
    samples = pd.read_csv(args.csv)

    variables = samples.columns.tolist()
    exclude = ["sampleNumber", "sampleName", "filePath", "genome", "controlSampleName"]
    [variables.pop(variables.index(exc)) for exc in exclude if exc in variables]

    cols = ["unpairedReadsExamined", "readPairsExamined", "unmappedReads", "unpairedReadDuplicates", "readPairDuplicates", "readPairOpticalDuplicates", "percentDuplication", "estimatedLibrarySize"]
    for col in cols:
        samples[col] = None

    for sample in xrange(len(samples)):
        # Sample name
        # if sampleName is not provided, use a concatenation of several variable (excluding longest)
        if str(samples["sampleName"][sample]) != "nan":
            sampleName = samples["sampleName"][sample]
        else:
            sampleName = string.join([str(samples[var][sample]) for var in variables], sep="_")
            logger.debug("No sample name provided, using concatenation of variables supplied")

        # Get duplicates
        dups = pd.read_csv(os.path.join(resultsDir, sampleName, ".duplicates.txt"), sep="\t", comment="#", header=1)
        dups.dropna(thresh=len(dups.columns) - 1, inplace=True)

        # Add values to sample sheet
        samples[[cols]] = dups.drop(["LIBRARY"], axis=1)

    # write annotation sheet with biological replicates
    samples.to_csv(os.path.join(projectDir, args.project_name + ".read_stats.csv"), index=False)

    logger.debug("Finished getting read statistics.")

    return (logger,
            os.path.join(os.getcwd(), args.project_name + ".log"),
            os.path.join(projectDir, "runs", args.project_name + ".log"))


def analyse(args, logger):
    logger.info("Starting sample analysis.")

    logger.debug("Checking project directories exist and creating if not.")
    htmlDir, projectDir, dataDir, resultsDir, urlRoot = checkProjectDirs(args, logger)

    # Paths to static files on the cluster

    # Other static info
    tagment = [
        "DNASE", "DNASESEQ", "DNASE-SEQ", "DHS", "DHS-SEQ", "DHSSEQ",
        "ATAC", "ATAC-SEQ", "ATACSEQ",
        "CM"
    ]
    histones = ["H2A", "H2B", "H3", "H4"]
    broadFactors = [
        "H3K27ME1", "H3K27ME2", "H3K27ME3",
        "H3K36ME1", "H3K36ME2", "H3K36ME3",
        "H3K9ME1", "H3K9ME2", "H3K9ME3",
        "H3K72ME1", "H3K72ME2", "H3K72ME3"
    ]

    tssFiles = {
        "hg19": "/fhgfs/groups/lab_bock/arendeiro/share/GRCh37_hg19_refSeq.tss.bed",
        "mm10": "/fhgfs/groups/lab_bock/arendeiro/share/GRCm38_mm10_refSeq.tss.bed",
        "dr7": "/fhgfs/groups/lab_bock/arendeiro/share/GRCh37_hg19_refSeq.tss.bed"
    }

    # Parse sample information
    args.csv = os.path.abspath(args.csv)

    # check if exists and is a file
    if not os.path.isfile(args.csv):
        logger.error("Sample annotation '%s' does not exist, or user has no read access." % args.csv)
        sys.exit(1)

    # read in
    samples = pd.read_csv(args.csv)

    # TODO:
    # Perform checks on the variables given
    # (e.g. columns existing, bams existing)

    # start pipeline
    projectName = string.join([args.project_name, time.strftime("%Y%m%d-%H%M%S")], sep="_")

    # Preprocess samples
    variables = samples.columns.tolist()
    exclude = ["sampleNumber", "sampleName", "filePath", "genome", "controlSampleName"]
    [variables.pop(variables.index(exc)) for exc in exclude if exc in variables]

    # track jobs to submit
    jobs = dict()

    for sample in xrange(len(samples)):
        # Sample name
        # if sampleName is not provided, use a concatenation of several variable (excluding longest)
        if str(samples["sampleName"][sample]) != "nan":
            sampleName = samples["sampleName"][sample]
        else:
            sampleName = string.join([str(samples[var][sample]) for var in variables], sep="_")
            logger.debug("No sample name provided, using concatenation of variables supplied")

        # Control name
        control = False
        # get file path of sample which this sample's controlName is pointing to
        ctrlField = samples[samples['sampleName'].isin([samples['controlSampleName'][sample]])]['filePath']

        # if there is only one record, use that as control
        if len(ctrlField) == 1:
            control = True
            controlName = samples.ix[sample]["controlSampleName"]
            controlBam = os.path.abspath(ctrlField.values[0])

        # check if sample is tagmented or not:
        tagmented = True if samples["technique"][sample] in tagment else False

        # Is it a histone?
        histone = True if any([i in samples["ip"][sample] for i in histones]) else False

        jobName = projectName + "_" + sampleName

        tempFiles = list()

        # assemble commands
        jobCode = slurmHeader(
            jobName=jobName,
            output=os.path.join(projectDir, "runs", jobName + ".slurm.log"),
            queue=args.queue,
            time=args.time,
            cpusPerTask=args.cpus,
            memPerCpu=args.mem,
            userMail=args.user_mail
        )
        if args.stage in ["all", "callpeaks"] and control:
            if args.peak_caller == "macs2":
                # make dir for output
                if not os.path.exists(os.path.join(dataDir, "peaks", sampleName)):
                    os.makedirs(os.path.join(dataDir, "peaks", sampleName))

                # For point-source factors use default settings
                # For broad factors use broad settings
                jobCode += macs2CallPeaks(
                    treatmentBam=os.path.abspath(samples['filePath'][sample]),
                    controlBam=controlBam,
                    outputDir=os.path.join(dataDir, "peaks", sampleName),
                    sampleName=sampleName,
                    genome=samples['genome'][sample],
                    broad=False if samples["ip"][sample].upper() not in broadFactors else True
                )
            elif args.peak_caller == "spp":
                # For point-source factors use default settings
                # For broad factors use broad settings
                jobCode += sppCallPeaks(
                    treatmentBam=os.path.abspath(samples['filePath'][sample]),
                    controlBam=controlBam,
                    treatmentName=sampleName,
                    controlName=controlName,
                    outputDir=os.path.join(dataDir, "peaks", sampleName),
                    broad=False if samples["ip"][sample].upper() not in broadFactors else True,
                    cpus=args.cpus
                )

        if args.stage in ["all", "findmotifs"] and control:
            # make dir for motifs
            if not os.path.exists(os.path.join(dataDir, "motifs", sampleName)):
                    os.makedirs(os.path.join(dataDir, "motifs", sampleName))

            if not histone:
                # For TFs, find the "self" motif
                jobCode += homerFindMotifs(
                    peakFile=os.path.join(dataDir, "peaks", sampleName, sampleName + "_peaks.narrowPeak"),
                    genome=samples["genome"][sample],
                    outputDir=os.path.join(dataDir, "motifs", sampleName),
                    size="50",
                    length="8,10,12,14,16",
                    n_motifs=8
                )
                # For TFs, find co-binding motifs (broader region)
                jobCode += homerFindMotifs(
                    peakFile=os.path.join(dataDir, "peaks", sampleName, sampleName + "_peaks.narrowPeak"),
                    genome=samples["genome"][sample],
                    outputDir=os.path.join(dataDir, "motifs", sampleName + "_cobinders"),
                    size="200",
                    length="8,10,12,14,16",
                    n_motifs=12
                )
            else:
                # For histones, use a broad region to find motifs
                jobCode += homerFindMotifs(
                    peakFile=os.path.join(dataDir, "peaks", sampleName, sampleName + "_peaks.narrowPeak"),
                    genome=samples["genome"][sample],
                    outputDir=os.path.join(dataDir, "motifs", sampleName),
                    size="1000",
                    length="8,10,12,14,16",
                    n_motifs=20
                )

        if args.stage in ["all", "centerpeaks"] and control:
            # TODO:
            # right now this assumes peaks were called with MACS2
            # figure a way of magetting the peak files withough using the peak_caller option
            # for that would imply taht option would be required when selecting this stage
            jobCode += centerPeaksOnMotifs(
                peakFile=os.path.join(dataDir, "peaks", sampleName, sampleName + "_peaks.narrowPeak"),
                genome=samples["genome"][sample],
                windowWidth=args.peak_window_width,
                motifFile=os.path.join(dataDir, "motifs", sampleName, "homerResults", "motif1.motif"),
                outputBed=os.path.join(dataDir, "peaks", sampleName, sampleName + "_peaks.motifCentered.bed")
            )

        if args.stage in ["all", "annotatepeaks"] and control:
            # TODO:
            # right now this assumes peaks were called with MACS2
            # figure a way of getting the peak files withough using the peak_caller option
            # for that would imply taht option would be required when selecting this stage
            jobCode += AnnotatePeaks(
                peakFile=os.path.join(dataDir, "peaks", sampleName, sampleName + "_peaks.narrowPeak"),
                genome=samples["genome"][sample],
                motifFile=os.path.join(dataDir, "motifs", sampleName, "homerResults", "motif1.motif"),
                outputBed=os.path.join(dataDir, "peaks", sampleName, sampleName + "_peaks.motifAnnotated.bed")
            )

        if args.stage in ["all", "peakanalysis"] and control:
            jobCode += peakAnalysis(
                inputBam=os.path.abspath(samples.ix[sample]['filePath']),
                peakFile=os.path.join(dataDir, "peaks", sampleName, sampleName + "_peaks.motifCentered.bed"),
                plotsDir=os.path.join(resultsDir, 'plots'),
                windowWidth=2000,
                fragmentsize=1 if tagmented else 50,  # change this to actual read length
                genome=samples.ix[sample]['genome'],
                n_clusters=5,
                strand_specific=True,
                duplicates=True
            )

        if args.stage in ["all", "tssanalysis"] and control:
            jobCode += tssAnalysis(
                inputBam=os.path.abspath(samples.ix[sample]['filePath']),
                tssFile=tssFiles[samples.ix[sample]['genome']],
                plotsDir=os.path.join(resultsDir, 'plots'),
                windowWidth=2000,
                fragmentsize=1 if tagmented else 50,  # change this to actual read length
                genome=samples.ix[sample]['genome'],
                n_clusters=5,
                strand_specific=True,
                duplicates=True
            )

        # if args.stage in ["all", "footprints"] and control and tagmented:
        #     jobCode += footprintAnalysis(
        #         os.path.join(dataDir, "peaks", sampleName, sampleName + "_peaks.motifCentered.bed"),
        #         os.path.join(dataDir, "peaks", sampleName, sampleName + "_peaks.motifAnnotated.bed")
        #     )

        # finish job
        jobCode += slurmFooter()
        jobs[jobName] = jobCode

    # Submit jobs to slurm
    for jobName, jobCode in jobs.items():
        # Output file name
        jobFile = os.path.join(projectDir, "runs", jobName + ".sh")

        with open(jobFile, 'w') as handle:
            handle.write(textwrap.dedent(jobCode))

        # Submit to slurm
        if not args.dry_run:
            logger.info("Submitting jobs to slurm")
            status = slurmSubmitJob(jobFile)

            if status != 0:
                logger.error("Slurm job '%s' not successfull" % jobFile)
                sys.exit(1)
            logger.debug("Project '%s'submission finished successfully." % args.project_name)

    logger.debug("Finished comparison")

    return (logger,
            os.path.join(os.getcwd(), args.project_name + ".log"),
            os.path.join(projectDir, "runs", args.project_name + ".log"))


def compare(args, logger):
    logger.info("Starting sample comparison.")

    logger.debug("Checking project directories exist and creating if not.")
    htmlDir, projectDir, dataDir, resultsDir, urlRoot = checkProjectDirs(args, logger)

    # Paths to static files on the cluster

    # Other static info

    # Parse sample information
    args.csv = os.path.abspath(args.csv)

    # check if exists and is a file
    if not os.path.isfile(args.csv):
        logger.error("Sample annotation '%s' does not exist, or user has no read access." % args.csv)
        sys.exit(1)

    # read in
    samples = pd.read_csv(args.csv)

    # TODO:
    # Perform checks on the variables given
    # (e.g. columns existing, bams existing)

    # start pipeline
    projectName = string.join([args.project_name, time.strftime("%Y%m%d-%H%M%S")], sep="_")

    variables = samples.columns.tolist()
    exclude = ["sampleNumber", "sampleName", "filePath", "genome", "controlSampleName"]
    [variables.pop(variables.index(exc)) for exc in exclude if exc in variables]

    # track jobs to submit
    jobs = dict()

    # diffBind
    if args.stage in ["all", "diffbind"]:

        # Submit separate jobs for each genome
        genome_groups = samples.groupby(["genome"]).groups

        for genome in genome_groups.keys():

            # Separate comparison per IPed factor
            df = samples.ix[genome_groups[genome]]
            IP_groups = df.groupby(["ip"]).groups

            for IP in IP_groups.keys():
                if IP.upper() in ["INPUT", "IGG"]:  # skip groups with control
                    continue

                jobName = projectName + "_" + "diffBind_{0}_{1}".format(genome, IP)
                diffBindSheetFile = os.path.join(projectDir, "runs", jobName + ".csv")

                # make diffBind csv file, save it
                empty = makeDiffBindSheet(samples, samples.ix[IP_groups[IP]], os.path.join(dataDir, "peaks"), diffBindSheetFile)

                if not empty:
                    # create job
                    jobCode = slurmHeader(
                        jobName=jobName,
                        output=os.path.join(projectDir, "runs", jobName + ".slurm.log"),
                        queue=args.queue,
                        time=args.time,
                        cpusPerTask=args.cpus,
                        memPerCpu=args.mem,
                        userMail=args.user_mail
                    )
                    jobCode += diffBind(
                        inputCSV=diffBindSheetFile,
                        jobName=jobName,
                        plotsDir=os.path.join(resultsDir, 'plots')
                    )
                    jobCode += slurmFooter()
                    jobs[jobName] = jobCode

    # Submit job for all samples together
    jobName = projectName + "_" + "sample_correlations"

    if args.stage in ["all", "correlations"]:
        # TODO:
        # Submit separate jobs for each genome
        jobCode = slurmHeader(
            jobName=jobName,
            output=os.path.join(projectDir, "runs", jobName + ".slurm.log"),
            queue=args.queue,
            time=args.time,
            cpusPerTask=args.cpus,
            memPerCpu=args.mem,
            userMail=args.user_mail
        )
        jobCode += plotCorrelations(
            inputBams=list(samples['filePath']),
            plotsDir=os.path.join(resultsDir, 'plots'),
            duplicates=args.duplicates,
            windowWidth=args.genome_window_width,
            fragmentSize=50,
            genome="hg19"
        )
        jobCode += slurmFooter()
        jobs[jobName] = jobCode

    # Submit jobs to slurm
    for jobName, jobCode in jobs.items():
        # Output file name
        jobFile = os.path.join(projectDir, "runs", jobName + ".sh")

        with open(jobFile, 'w') as handle:
            handle.write(textwrap.dedent(jobCode))

        # Submit to slurm
        if not args.dry_run:
            logger.info("Submitting jobs to slurm")
            status = slurmSubmitJob(jobFile)

            if status != 0:
                logger.error("Slurm job '%s' not successfull" % jobFile)
                sys.exit(1)
            logger.debug("Project '%s'submission finished successfully." % args.project_name)

    logger.debug("Finished comparison")

    return (logger,
            os.path.join(os.getcwd(), args.project_name + ".log"),
            os.path.join(projectDir, "runs", args.project_name + ".log"))


def makeDiffBindSheet(samples, df, peaksDir, sheetFile):
    """
    Prepares and saves a diffBind annotation sheet from a pandas.DataFrame annotation with standard columns.
    """
    for col in ["numberCells", "technique", "treatment", "patient"]:
        df[col] = df[col].apply(str)
    df["condition"] = df["numberCells"] + "_" + df["technique"]
    df["treatment"] = df["treatment"] + "_" + df["patient"]
    df["ControlID"] = list(samples["sampleName"][list(df["controlSampleName"] - 1)])
    df["bamControl"] = list(samples["filePath"][list(df["controlSampleName"] - 1)])
    df["Peaks"] = [os.path.join(peaksDir, sampleName, sampleName + "_peaks.narrowPeak") for sampleName in df["sampleName"]]
    df["PeakCaller"] = "macs"
    df = df[["sampleName", "cellLine", "ip", "condition", "treatment",
             "biologicalReplicate", "filePath", "ControlID", "bamControl", "Peaks", "PeakCaller"]]
    df.columns = ["SampleID", "Tissue", "Factor", "Condition", "Treatment",
                  "Replicate", "bamReads", "ControlID", "bamControl", "Peaks", "PeakCaller"]
    # exclude samples without matched control
    df = df[df["ControlID"].notnull()]
    # if not empty
    if not df.empty:
        # save as csv
        df.to_csv(sheetFile, index=False)  # check if format complies with DiffBind

    return df.empty


def checkSamplePE(sampleFile, n=1000):
    p = subprocess.Popen(['samtools', 'view', sampleFile], stdout=subprocess.PIPE)

    # Count paired alignments
    paired = 0
    while n > 0:
        if 1 & int(p.stdout.next().split("\t")[1]):  # check decimal flag contains 1 (paired)
            paired += 1
        n -= 1
    # If at least half is paired, return True
    if paired < (n / 2):
        return False
    else:
        return True


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
    {1}
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
    # Merging bam files from replicates
    echo "Merging bam files from replicates"

    java -Xmx4g -jar /cm/shared/apps/picard-tools/1.118/SamToFastq.jar \\
    INPUT={0} """.format(inputBam)
    if outputFastq2 is None and unpairedFastq is None:
        command += "FASTQ={0}".format(outputFastq)
    else:
        command += """FASTQ={0} \\
        SECOND_END_FASTQ={1} \\
        UNPAIRED_FASTQ={2}
        """.format(outputFastq, outputFastq2, unpairedFastq)

    return command


def trimAdaptersSE(inputFastq, outputFastq, cpus, adapters, log):
    command = """
    # Trimming adapters from sample
    echo "Trimming adapters from sample"
    module load trimmomatic/0.32

    java -Xmx4g -jar `which trimmomatic-0.32.jar` SE \\
    -threads {0} \\
    -trimlog {1} \\
    {2} \\
    {3} \\
    ILLUMINACLIP:{4}:1:40:15 \\
    LEADING:3 TRAILING:3 \\
    SLIDINGWINDOW:4:10 \\
    MINLEN:36

    """.format(cpus, log, inputFastq, outputFastq, adapters)

    return command


def trimAdaptersPE(inputFastq1, inputFastq2,
                   outputFastq1, outputFastq1unpaired,
                   outputFastq2, outputFastq2unpaired,
                   cpus, adapters, log):
    command = """
    # Trimming adapters from sample
    echo "Trimming adapters from sample"
    module load trimmomatic/0.32

    java -Xmx4g -jar `which trimmomatic-0.32.jar` PE \\
    -threads {0} \\
    -trimlog {1} \\
    {2} \\
    {3} \\
    {4} \\
    {5} \\
    {6} \\
    {7} \\
    ILLUMINACLIP:{8}:1:40:15 \\
    LEADING:3 TRAILING:3 \\
    SLIDINGWINDOW:4:10 \\
    MINLEN:36

    """.format(cpus, log,
               inputFastq1, inputFastq2,
               outputFastq1, outputFastq1unpaired,
               outputFastq2, outputFastq2unpaired,
               adapters)

    return command


def bowtie2Map(inputFastq, outputBam, genomeIndex, cpus, inputFastq2=None):
    outputBam = re.sub("\.bam$", "", outputBam)
    # Will only admit 500bp-long fragments. Change with --maxins option
    command = """
    # Mapping reads with Bowtie2
    echo "Mapping reads with Bowtie2"
    module load bowtie/2.2.3
    module load samtools

    bowtie2 --very-sensitive -p {0} \\
    -x {1} """.format(cpus, genomeIndex)
    if inputFastq2 is None:
        command += "{0}".format(inputFastq)
    else:
        command += " -1 {0} -2 {1}".format(inputFastq, inputFastq2)
    command += """ | \\
    samtools view -S -b - | \\
    samtools sort - {3}

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
    samtools sort {1} {4}

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
    sambamba markdup -t {2} -r {0} {1}

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
    genomeCoverageBed -i stdin -bg -g {1} > {2}.cov

    bedGraphToBigWig {2}.cov {1} {3}

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
    genomeCoverageBed -i stdin -bg -g {1} > {2}.cov

    bedGraphToBigWig {2}.cov {1} {3}

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
    echo 'track type=bigWig name="{0} {1}" description="{0} {1}" """.format(sampleName, fivePrime)
    command += """height=32 visibility=full maxHeightPixels=32:32:25 bigDataUrl={0} color={1}' >> {2}

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


def macs2CallPeaks(treatmentBam, controlBam, outputDir, sampleName, genome, broad=False):
    if genome == "hg19":
        genome = "hs"
    elif genome == "mm10":
        genome = "mm"
    elif genome == "dr7":
        raise ValueError("Genome dr7 not yet supported for peak calling with MACS2.")

    if not broad:
        # Check wth Christian if fragment length should be different for ChIP or CM
        command = """
    # Call peaks with MACS2
    echo "Calling peaks with MACS2"

    macs2 callpeak \\
    -t {0} \\
    -c {1} \\
    --bw 200 \\
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


def plotCorrelations(inputBams, plotsDir, duplicates=False,
                     windowWidth=1000, fragmentSize=50, genome="hg19"):
    command = """
    # Plot correlations
    echo "plotting correlations"

    source ~/venv/bin/activate

    python {5}/lib/correlations.py {0} \\
    {1} \\
    --window-width {2} \\
    --fragment-size {3} \\
    --genome {4}

    """.format(plotsDir, " ".join(["%s"] * len(inputBams)) % tuple(inputBams),
               windowWidth, fragmentSize, genome,
               os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
               )

    if duplicates:
        command += "--duplicates"

    command += "\n    deactivate\n"

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


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
