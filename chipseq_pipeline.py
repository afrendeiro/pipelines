#!/usr/bin/env python

"""
ChIP-seq pipeline

    ## BUGS
    Keep track of original bam file location (1 technical replicate)
    Info "exiting" after copying log
    chmod 755 parent bigWig folder

    UCSC link for mouse is: http://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&hgt.customText=

    UCSC track description 5prime add 5prime

    Take care of additional requirements: UCSCutils, ...


    ## FURTHER TO IMPLEMENT
    check annotation sheets are right
    call footprints
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
                                               "indexbam", "qc", "maketracks", "mergereplicates"],
                                      help='Run only these stages. Default=all.', type=str)
    # comparison
    comparison_subparser = subparser.add_parser("comparison")
    comparison_subparser.add_argument(dest='project_name', help="Project name.", type=str)
    comparison_subparser.add_argument(dest='csv', help='CSV file with sample annotation.', type=str)
    comparison_subparser.add_argument('-s', '--stage', default="all", dest='stage',
                                      choices=["all", "callpeaks", "findmotifs", "centerpeaks",
                                               "annotatepeaks", "peakanalysis", "footprints", "correlations"],
                                      help='Run only these stages. Default=all.', type=str)
    comparison_subparser.add_argument('--peak-caller', default="macs2", choices=["macs2", "spp"],
                                      dest='peak_caller', help='Peak caller to use. Default=macs2.', type=str)
    comparison_subparser.add_argument('--peak-window-width', default=2000,
                                      dest='peak_window_width', help='Width of window around peak motifs. Default=2000.', type=int)
    comparison_subparser.add_argument('--duplicates', dest='duplicates', action='store_true',
                                      help='Allow duplicates in coorelation analysis. Default=False.')
    comparison_subparser.add_argument('--genome-window-width', default=1000,
                                      dest='genome_window_width', help='Width of window to make genome-wide correlations. Default=1000.', type=int)

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
        preprocess(args, logger)
    elif args.command == "comparison":
        comparison(args, logger)

    # Exit
    logger.info("Finished and exiting.")
    sys.exit(0)


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
        os.path.join(dataDir, "peaks"),
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

    return (htmlDir, projectDir, dataDir, resultsDir, urlRoot)


def checkTechnicalReplicates(samples):
    """
    returns dictionary with entries and list of pd.Series with sample info
    samples - a pandas.DataFrame with sample info.

    """
    # Check if there are technical replicates
    variables = samples.columns.tolist()
    exclude = ["sampleNumber", "sampleName", "technicalReplicate", "filePath", "controlSampleNumber"]
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
        "hg19": os.path.join(genomeFolder, "hg19/forBowtie2/hg19"),
        "mm10": os.path.join(genomeFolder, "mm10/forBowtie2/mm10"),
        "dr7": os.path.join(genomeFolder, "dr7/forBowtie2/dr7")
    }
    genomeSizes = {
        "hg19": os.path.join(genomeFolder, "hg19/hg19_chromlengths.txt"),
        "mm10": "/home/arendeiro/mm10.chrom.sizes",
        "dr7": "/home/arendeiro/danRer7.chrom.sizes"
    }
    adapterFasta = "/fhgfs/groups/lab_bock/shared/chipmentation.fa"

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
        exclude = ["sampleNumber", "sampleName", "filePath", "genome", "tagmented", "controlSampleNumber"]
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

        # check if sample is tagmented or not:
        # todo: get tagmented from technique
        tagmented = biologicalReplicates[sample][0]["tagmented"] == "yes" or biologicalReplicates[sample][0]["tagmented"] == 1

        if not tagmented:
            bam = os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2")
        else:
            bam = os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.shifted")

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
                    outputBam=os.path.join(dataDir, "raw", sampleName + ".bam")
                )

            # convert bam to fastq
            if args.stage in ["all", "bam2fastq"]:
                if len(biologicalReplicates.values()[sample]) == 1:
                    jobCode += bam2fastq(
                        inputBam=biologicalReplicates.values()[sample][0]["filePath"],
                        outputFastq=os.path.join(dataDir, "fastq", sampleName + ".fastq")
                    )
                elif len(biologicalReplicates.values()[sample]) >= 1:
                    jobCode += bam2fastq(
                        inputBam=os.path.join(dataDir, "raw", sampleName + ".bam"),
                        outputFastq=os.path.join(dataDir, "fastq", sampleName + ".fastq")
                    )
                tempFiles.append(os.path.join(dataDir, "fastq", sampleName + ".fastq"))
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
                outputFastq=os.path.join(dataDir, "fastq", sampleName + ".trimmed.fastq"),
                adapters=adapterFasta
            )
            tempFiles.append(os.path.join(dataDir, "fastq", sampleName + ".trimmed.fastq"))
        if args.stage in ["all", "mapping"]:
            if samples["genome"][sample] not in genomeIndexes:
                logger.error("Sample %s has unsuported genome index: %s" % (sampleName, samples["genome"][sample]))
                sys.exit(1)
            jobCode += bowtie2Map(
                inputFastq=os.path.join(dataDir, "fastq", sampleName + ".trimmed.fastq"),
                outputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.bam"),
                genomeIndex=genomeIndexes[samples["genome"][sample]],
                cpus=args.cpus
            )
            tempFiles.append(os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.bam"))
        if args.stage in ["all", "shiftreads"]:
            if tagmented:
                # TODO:
                # Get correct relative path.
                # Relative to location of *this* file rather than relative to cwd.
                jobCode += shiftReads(
                    inputBam=os.path.join(dataDir, "mapped", sampleName + ".trimmed.bowtie2.bam"),
                    outputBam=bam + ".bam"
                )
        if args.stage in ["all", "markduplicates"]:
            jobCode += markDuplicates(
                inputBam=bam + ".bam",
                outputBam=bam + ".dups.bam",
                metricsFile=os.path.join(resultsDir, sampleName + ".duplicates.txt")#,
                #tempDir=
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
                tagmented=False
            )
            jobCode += addTrackToHub(
                sampleName=sampleName,
                trackURL=urlRoot + sampleName + ".bigWig",
                trackHub=os.path.join(htmlDir, "trackHub.txt")
            )
            if tagmented:
                jobCode += bamToUCSC(
                    inputBam=bam + ".dups.bam",
                    outputBigWig=os.path.join(htmlDir, sampleName + ".5prime.bigWig"),
                    genomeSizes=genomeSizes[samples["genome"][sample]],
                    tagmented=True
                )
                jobCode += addTrackToHub(
                    sampleName=sampleName,
                    trackURL=urlRoot + sampleName + ".5prime.bigWig",
                    trackHub=os.path.join(htmlDir, "trackHub.txt"),
                    fivePrime="5prime"
                )
            # TODO: separate this per genome
            linkToTrackHub(
                trackHubURL="{0}/{1}/bigWig/trackHub.txt".format(args.url_root, args.project_name),
                fileName=os.path.join(projectDir, "ucsc_tracks.html"),
                genome='human'
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
    df["controlSampleNumber"] = None
    df.to_csv(os.path.join(projectDir, args.project_name + ".biol_replicates.annotation_sheet.csv"), index=False)

    # Copy log to projectDir
    shutil.copy2(
        os.path.join(os.getcwd(), args.project_name + ".log"),
        os.path.join(projectDir, "runs", args.project_name + ".log")
    )
    logger.debug("Copied log file to project directory '%s'" % os.path.join(projectDir, "runs"))


def comparison(args, logger):
    logger.info("Starting sample comparison.")

    logger.debug("Checking project directories exist and creating if not.")
    htmlDir, projectDir, dataDir, resultsDir, urlRoot = checkProjectDirs(args, logger)

    # Paths to static files on the cluster

    # Other static info
    broadFactors = ["H3K27me3", "H3K36me3", "H3K9me3"]

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
    exclude = ["sampleNumber", "sampleName", "filePath", "genome", "tagmented", "controlSampleNumber"]
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
        if not np.isnan(samples.ix[sample]["controlSampleNumber"]):
            control = True
            controlIdx = int(samples.ix[sample]["controlSampleNumber"]) - 1

            # if sampleName is not provided, use a concatenation of several variable (excluding longest)
            if str(samples["sampleName"][controlIdx]) != "nan":
                controlName = samples["sampleName"][controlIdx]
            else:
                controlName = string.join([str(samples[var][controlIdx]) for var in variables], sep="_")
                logger.debug("No sample name provided, using concatenation of variables supplied")

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

                if samples["ip"][sample] not in broadFactors:
                    jobCode += macs2CallPeaks(
                        treatmentBam=os.path.abspath(samples.ix[sample]['filePath']),
                        controlBam=os.path.abspath(samples.ix[controlIdx]['filePath']),
                        outputDir=os.path.join(dataDir, "peaks", sampleName),
                        sampleName=sampleName,
                        broad=False
                    )
                else:
                    jobCode += macs2CallPeaks(treatmentBam=os.path.abspath(samples.ix[sample]['filePath']),
                                              controlBam=os.path.abspath(samples.ix[controlIdx]['filePath']),
                                              outputDir=os.path.join(dataDir, "peaks", sampleName),
                                              sampleName=sampleName,
                                              broad=True
                                              )
            elif args.peak_caller == "spp":
                jobCode += sppCallPeaks(
                    treatmentBam=os.path.abspath(samples.ix[sample]['filePath']),
                    controlBam=os.path.abspath(samples.ix[controlIdx]['filePath']),
                    treatmentName=sampleName,
                    controlName=controlName,
                    outputDir=os.path.join(dataDir, "peaks", sampleName),
                    cpus=args.cpus
                )

        if args.stage in ["all", "findmotifs"] and control:
            if not os.path.exists(os.path.join(dataDir, "motifs", sampleName)):
                    os.makedirs(os.path.join(dataDir, "motifs", sampleName))

            jobCode += homerFindMotifs(peakFile=os.path.join(dataDir, "peaks", sampleName, sampleName + "_peaks.narrowPeak"),
                                       genome=samples["genome"][sample],
                                       outputDir=os.path.join(dataDir, "motifs", sampleName),
                                       size="150",
                                       length="8,10,12"
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
            # figure a way of magetting the peak files withough using the peak_caller option
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
                fragmentsize=1 if samples.ix[sample]['tagmented'] else 50,
                genome=samples.ix[sample]['genome'],
                n_clusters=5,
                strand_specific=True,
                duplicates=True
            )

        # if args.stage in ["all", "footprints"] and control:
        #     jobCode += footprintAnalysis(
        #         os.path.join(dataDir, "peaks", sampleName, sampleName + "_peaks.motifCentered.bed"),
        #         os.path.join(dataDir, "peaks", sampleName, sampleName + "_peaks.motifAnnotated.bed")
        #     )

        # finish job
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

    # Copy log to projectDir
    shutil.copy2(
        os.path.join(os.getcwd(), args.project_name + ".log"),
        os.path.join(projectDir, "runs", args.project_name + ".log")
    )
    logger.debug("Copied log file to project directory '%s'" % os.path.join(projectDir, "runs"))


def slurmHeader(jobName, output, queue="shortq", ntasks=1, time="10:00:00", cpusPerTask=16, memPerCpu=2000, nodes=1, userMail=""):
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

    """.format(queue, ntasks, time, cpusPerTask, memPerCpu, nodes, jobName, output, userMail)

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
    module load samtools

    samtools merge {0} {1}

    """.format(outputBam, (" ".join(["%s"] * len(inputBams))) % tuple(inputBams))

    return command


def bam2fastq(inputBam, outputFastq):
    command = """
    # Converting original Bam file to Fastq
    echo "Converting original Bam file to Fastq"
    module load bamtools

    bamtools convert -in {0} -format fastq -out {1}

    """.format(inputBam, outputFastq)

    return command


def fastqc(inputFastq, outputDir):
    command = """
    # Measuring sample quality with Fastqc
    echo "Measuring sample quality with Fastqc"
    module load java/jdk/1.7.0_65
    module load FastQC/0.11.2

    fastqc --noextract --outdir {0} {1}

    """.format(outputDir, inputFastq)

    return command


def trimAdapters(inputFastq, outputFastq, adapters):
    command = """
    # Trimming adapters from sample
    echo "Trimming adapters from sample"
    module load trimmomatic/0.32

    java -jar `which trimmomatic-0.32.jar` SE {0} {1} \\
    ILLUMINACLIP:{2}:1:40:15 \\
    LEADING:3 TRAILING:3 \\
    SLIDINGWINDOW:4:10 \\
    MINLEN:36

    """.format(inputFastq, outputFastq, adapters)

    return command


def bowtie2Map(inputFastq, outputBam, genomeIndex, cpus):
    outputBam = re.sub("\.bam$", "", outputBam)
    command = """
    # Mapping reads with Bowtie2
    echo "Mapping reads with Bowtie2"
    module load bowtie/2.2.3
    module load samtools

    bowtie2 --very-sensitive -p {0} -x {1} {2} | \\
    samtools view -S -b - | \\
    samtools sort - {3}

    """.format(cpus, genomeIndex, inputFastq, outputBam)

    return command


def shiftReads(inputBam, outputBam):
    outputBam = re.sub("\.bam$", "", outputBam)
    # TODO:
    # Implement read shifting with HTSeq or Cython
    command = """
    # Shifting read of tagmented sample
    echo "Shifting read of tagmented sample"
    module load samtools
    module load python

    samtools view -h {0} | \\
    python {1}/lib/shift_reads.py | \\
    samtools view -S -b - | \\
    samtools sort - {2}

    """.format(inputBam, os.path.abspath(os.path.dirname(os.path.realpath(__file__))), outputBam)

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


def bamToUCSC(inputBam, outputBigWig, genomeSizes, tagmented=False):
    transientFile = os.path.abspath(re.sub("\.bigWig", "", outputBigWig))
    if not tagmented:
        command = """
    # Making bigWig tracks from bam file
    echo "making bigWig tracks from bam file"
    module load bedtools

    bedtools bamtobed -i {0} | \\
    genomeCoverageBed -i stdin -bg -g {1} > {2}.cov

    bedGraphToBigWig {2}.cov {1} {3}

    # remove cov file
    if [[ -s {2}.cov ]]
        then
        rm {2}.cov
    fi

    chmod 755 {3}

    """.format(inputBam, genomeSizes, transientFile, outputBigWig)
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

    """.format(inputBam, genomeSizes, transientFile, outputBigWig, os.path.abspath(os.path.dirname(os.path.realpath(__file__))))

    return command


def addTrackToHub(sampleName, trackURL, trackHub, fivePrime=""):
    command = """
    # Adding track to TrackHub
    echo "Adding track to TrackHub"
    echo 'track type=bigWig name="{0} {3}" description="{0} {3}" visibility=3 bigDataUrl={1}' >> {2}

    chmod 755 {2}

    """.format(sampleName, trackURL, trackHub, fivePrime)

    return command


def linkToTrackHub(trackHubURL, fileName, genome):
    html = """
    <html>
        <head>
            <meta http-equiv="refresh" content="0; url=http://genome.ucsc.edu/cgi-bin/hgTracks?org={genome}&hgt.customText={trackHubURL}" />
        </head>
    </html>
    """.format(trackHubURL=trackHubURL, genome=genome)

    with open(fileName, 'w') as handle:
        handle.write(textwrap.dedent(html))


def macs2CallPeaks(treatmentBam, controlBam, outputDir, sampleName, broad=False):
    if not broad:
        command = """
    # Calling peaks with MACS2
    echo "calling peaks with MACS2"

    macs2 callpeak -t {0} -c {1} \\
    --bw 200 \\
    -g hs -n {2} --outdir {3}

    """.format(treatmentBam, controlBam, sampleName, outputDir)
    else:
        command = """
    # Call peaks with MACS2

    macs2 callpeak -B -t {0} -c {1} \\
    --broad --nomodel --extsize 200 --pvalue 1e-3 \\
    -g hs -n {2} --outdir {3}

    """.format(treatmentBam, controlBam, sampleName, outputDir)

    return command


def sppCallPeaks(treatmentBam, controlBam, treatmentName, controlName, outputDir, cpus):
    command = """
    # Calling peaks with SPP
    echo "calling peaks with SPP"
    module load R

    Rscript {6}/lib/spp_peak_calling.R {0} {1} {2} {3} {4} {5}

    """.format(treatmentBam, controlBam, treatmentName, controlName, cpus,
               os.path.abspath(os.path.dirname(os.path.realpath(__file__))))

    return command


def homerFindMotifs(peakFile, genome, outputDir, size=150, length="8,10,12"):
    command = """
    # Find motifs with Homer
    echo "finding motifs with Homer"

    findMotifsGenome.pl {0} {1} {2} -mask -size {3} -len {4}

    """.format(peakFile, genome, outputDir, size, length)

    return command


def AnnotatePeaks(peakFile, genome, motifFile, outputBed):
    command = """
    # Annotate peaks with motif score
    echo "annotating peaks with motif score"

    annotatePeaks.pl {0} {1} -mask -mscore -m {2} | \\
    tail -n +2 | cut -f 1,5,22  > {3}
    """.format(peakFile, genome, motifFile, outputBed)

    return command


def centerPeaksOnMotifs(peakFile, genome, windowWidth, motifFile, outputBed):
    command = """
    # Center peaks on motif
    echo "centering peaks on motif"

    annotatePeaks.pl {0} {1} -size {2} -center {3} | \\
    awk -v OFS='\\t' '{{print $2, $3, $4, $1, $6, $5}}' | \\
    python {4}/lib/fix_bedfile_genome_boundaries.py | \\
    sortBed > {5}
    """.format(peakFile, genome, windowWidth, motifFile, os.path.abspath(os.path.dirname(os.path.realpath(__file__))), outputBed)

    return command


def peakAnalysis(inputBam, peakFile, plotsDir, windowWidth, fragmentsize, genome, n_clusters, strand_specific, duplicates):
    command = """
    # Analyse peak profiles
    echo "Analysing peak profiles"

    {0}/lib/peaks_analysis.py {1} {2} {3} --window-width {4} --fragment-size {5} --genome {6} --n_clusters {7} """.format(
        os.path.abspath(os.path.dirname(os.path.realpath(__file__))),
        inputBam, peakFile, plotsDir, windowWidth, fragmentsize, genome, n_clusters
    )
    if strand_specific:
        command += "--strand-specific "
    if duplicates:
        command += "--duplicates "

    return command


def footprintAnalysis():
    raise NotImplementedError("Function not implemented yet.")


def plotCorrelations(inputBams, plotsDir, duplicates=False, windowWidth=1000, fragmentSize=50, genome="hg19"):
    command = """
    # Plot correlations
    echo "plotting correlations"

    source ~/venv/bin/activate

    python {5}/lib/correlations.py {0} {1} --window-width {2} --fragment-size {3} --genome {4}

    """.format(plotsDir, " ".join(["%s"] * len(inputBams)) % tuple(inputBams),
               windowWidth, fragmentSize, genome,
               os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
               )

    if duplicates:
        command += "--duplicates"

    command += "\n    deactivate\n"

    return command


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(0)
