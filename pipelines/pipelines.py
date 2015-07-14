#!/usr/bin/env python

"""
pipelines
=========

Project management and Sample loop.
"""

from argparse import ArgumentParser
from pipelines import Project
from pipelines import toolkit as tk
import cPickle as pickle
import os
import pandas as pd
import re
import subprocess
import sys
import textwrap
import time

__author__ = "Andre Rendeiro"
__copyright__ = "Copyright 2015, Andre Rendeiro"
__credits__ = []
__license__ = "GPL2"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"


def main():
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="pipelines",
        description="pipelines. Project management and sample loop."
    )
    parser = addArgs(parser)

    # Parse
    args = parser.parse_args()

    # Start project
    prj = Project(args.project_name)
    prj.addSampleSheet(args.csv)

    # Start main function
    if args.stats:
        readStats(args, prj)
    elif args.compare:
        compare(args, prj)
    else:
        sampleLoop(args, prj)

    # Exit
    print("Finished and exiting.")
    sys.exit(0)


def addArgs(parser):
    """
    Options for project and pipelines.
    """
    # Project
    parser.add_argument(dest="project_name", help="Project name.", type=str)
    parser.add_argument(dest="csv", help="CSV file with sample annotation.", type=str)  # improvement: check project dirs for csv

    # Behaviour
    parser.add_argument("--stats", dest="stats", action="store_true",
                        help="Do not run pipelines, but gather stats on produced files.")
    parser.add_argument("--compare", dest="compare", action="store_true",
                        help="Do not loop through samples, but perform comparisons betweem them.")
    parser.add_argument("-r", "--rm-tmp", dest="rm_tmp", action="store_true",
                        help="Remove intermediary files. If not it will preserve all intermediary files. Default=False.")
    parser.add_argument("--dry-run", dest="dry_run", action="store_true",
                        help="Dry run. Assemble commands, but do not submit jobs to slurm. Default=False.")
    parser.add_argument("--no-checks", dest="checks", action="store_false",
                        help="Don't check file existence and integrity. Default=False.")

    # Pypiper
    parser.add_argument("--overwrite", dest="recover", action="store_true",
                        help="Overwrite existing files. Default=False.")
    parser.add_argument("--fresh-start", dest="fresh", action="store_true",
                        help="Start from beginning of pipeline. Default=False.")
    parser.add_argument("--manual-clean", dest="manual_clean", action="store_true",
                        help="Manually clean temporary files. Default=False.")

    # Slurm-related
    parser.add_argument("-c", "--cpus", default=4, dest="cpus",
                        help="Number of CPUs to use. Default is specified in the pipeline config file.", type=int)
    parser.add_argument("-m", "--mem-per-cpu", default=4000, dest="mem",
                        help="Memory per CPU to use. Default is specified in the pipeline config file.", type=int)
    parser.add_argument("-q", "--queue", default="shortq", dest="queue",
                        choices=["develop", "shortq", "mediumq", "longq"],
                        help="Queue to submit jobs to. Default is specified in the pipeline config file.", type=str)
    parser.add_argument("-t", "--time", default="10:00:00", dest="time",
                        help="Maximum time for jobs to run. Default is specified in the pipeline config file.", type=str)
    parser.add_argument("-u", "--user-mail", default="mail@example.com", dest="user_mail",
                        help="User email.", type=str)

    # Preprocessing: trimming, mapping, etc...
    parser.add_argument("--trimmer", default="skewer", choices=["trimmomatic", "skewer"],
                        dest="trimmer", help="Trimmer to use. Default=skewer.", type=str)
    parser.add_argument("-i", "--max-insert-size", default=2000,
                        dest="maxinsert",
                        help="Maximum allowed insert size allowed for paired end mates. Default=2000.",
                        type=int)
    parser.add_argument("-Q", "--quality", default=30,
                        dest="quality",
                        help="Minimum read quality to keep. Default=30.",
                        type=int)

    # Further downstream
    parser.add_argument("--window-size", default=1000, dest="windowsize",
                        help="Window size used for genome-wide correlations. Default=1000.",
                        type=int)
    parser.add_argument("--peak-caller", default="macs2", choices=["macs2", "spp"],
                        dest="peak_caller", help="Peak caller to use. Default=macs2.", type=str)
    parser.add_argument("--peak-window-width", default=2000,
                        dest="peak_window_width",
                        help="Width of window around peak motifs. Default=2000.",
                        type=int)

    return parser


def sampleLoop(args, prj):
    """
    Loop through all samples and submit jobs to the pipeline under Slurm.

    :param args: Parsed ArgumentParser object.
    :type args: argparse.ArgumentParser
    :param prj: `Project` object.
    :type prj: pipelines.Project
    """

    print("Starting sample preprocessing into jobs.")

    # start pipeline
    runName = "_".join([prj.name, time.strftime("%Y%m%d-%H%M%S")])

    # add track headers to track hubs
    for genome in pd.Series([s.genome for s in prj.samples]).unique():
        with open(os.path.join(prj.dirs.html, "trackHub_{0}.txt".format(genome)), "w") as handle:
            handle.write("browser position {0}\n".format(prj.config["defaultposition"]))

    # Loop through samples, submit to corresponding job (preprocess, analyse)
    for sample in prj.samples:
        # get jobname
        jobName = "_".join([runName, sample.sampleName])

        # if unmappedBam is a list, add final "unmapped" attr to sample object
        if type(sample.unmappedBam) is list:
            sample.unmapped = os.path.join(sample.dirs.unmapped, sample.name + ".bam")

        # assemble command
        # slurm header
        jobCode = tk.slurmHeader(
            jobName=jobName,
            output=os.path.join(prj.dirs.logs, jobName + ".slurm.log"),
            queue=args.queue,
            time=args.time,
            cpusPerTask=args.cpus,
            memPerCpu=args.mem,
            userMail=args.user_mail
        )

        samplePickle = os.path.join(prj.dirs.pickles, jobName + ".pickle")
        # self reference the pickle file in its sample
        sample.pickle = samplePickle

        # If sample has control attribute, get that sample and pair them
        if hasattr(sample, "controlSampleName"):
            if type(sample.controlSampleName) == str:
                # Assign the sample with that name to ctrl
                ctrl = [s for s in prj.samples if s.sampleName == sample.controlSampleName]
                # if there is only one record, use that as control
                if len(ctrl) == 1:
                    sample.ctrl = ctrl[0]
                else:
                    # if not, process sample anyway, but without a matched control
                    print("Provided control sample name does not exist or is ambiguous: %s" % sample.controlSampleName)

        # save pickle with all objects (this time, 2nd element is a tuple!)
        pickle.dump((prj, sample, args), open(samplePickle, "wb"))

        # Actual call to pipeline
        technique = sample.technique.upper()
        if technique in prj.config["techniques"]["chipseq"]:
            jobCode += "chipseq-pipeline {0}\n".format(samplePickle)
        elif technique in prj.config["techniques"]["cm"]:
            jobCode += "chipseq-pipeline {0}\n".format(samplePickle)
        elif technique in prj.config["techniques"]["atacseq"]:
            jobCode += "atacseq-pipeline {0}\n".format(samplePickle)
        elif technique in prj.config["techniques"]["dnase"]:
            jobCode += "atacseq-pipeline {0}\n".format(samplePickle)
        elif technique in prj.config["techniques"]["quantseq"]:
            jobCode += "quantseq-pipeline {0}\n".format(samplePickle)
        else:
            raise TypeError("Sample is not in known sample class.")

        # Slurm footer
        jobCode += tk.slurmFooter()

        # Save code as executable
        jobFile = os.path.join(prj.dirs.executables, jobName + ".sh")
        with open(jobFile, 'w') as handle:
            handle.write(textwrap.dedent(jobCode))

        # Submit to slurm
        if not args.dry_run:
            status = tk.slurmSubmitJob(jobFile)

            if status != 0:
                print("Could not submit job '%s' to slurm." % jobFile)
                sys.exit(1)
            print("Submitted job to slurm: '%s'" % jobName)

        # Create link to trackHub in project folder
        tk.linkToTrackHub(
            trackHubURL=os.path.join(prj.dirs.html, "trackHub_{0}.txt".format(sample.genome)),
            fileName=os.path.join(prj.dirs.root, "ucsc_tracks_{0}.html".format(sample.genome)),
            genome=sample.genome
        )

    # write original annotation sheet to project folder
    # add field for manual sample-control pairing
    prj.sheet.df.controlSampleName = None
    prj.sheet.to_csv(os.path.join(prj.dirs.root, prj.name + ".annotation_sheet.csv"))

    print("Finished preprocessing")


def readStats(args, prj):
    """
    Given an annotation sheet with replicates, gets number of reads mapped, duplicates, etc...

    :param args: Parsed ArgumentParser object.
    :type args: argparse.ArgumentParser
    :param prj: `Project` object.
    :type prj: pipelines.Project
    """
    import re

    print("Starting sample read stats.")

    bowtieCols = ["readCount", "unpaired", "unaligned", "unique", "multiple", "alignmentRate"]

    samples = pd.DataFrame(index=["sampleName"] + bowtieCols + ["single-ends", "paired-ends", "duplicates", "NSC", "RSC", "qualityTag", "peakNumber", "FRiP"])

    for sample in prj.samples:
        sample = sample.asSeries()
        # Get alignment stats
        try:
            sample = sample.append(tk.parseBowtieStats(sample.alnRates))
        except:
            print("Record with alignment rates is empty or not found for sample %s" % sample.sampleName)
            pass

        # Get duplicate stats
        try:
            sample = sample.append(tk.parseDuplicateStats(sample.dupsMetrics))
        except:
            print("Record with duplicates is empty or not found for sample %s" % sample.sampleName)
            pass

        # Get NSC and RSC
        try:
            sample = sample.append(tk.parseQC(sample.sampleName, sample.qc))
        except:
            print("Record with quality control is empty or not found for sample %s" % sample.sampleName)
            pass

        # Count peak number (if peaks exist)
        if hasattr(sample, "peaks"):
            # and if sample has peaks
            if str(sample.peaks) != "nan":
                # if peak file exist
                if os.path.exists(sample.peaks):
                    sample = tk.getPeakNumber(sample)

        # Get FRiP from file (if exists) and add to sheet
        if hasattr(sample, "peaks"):
            # if sample has peaks
            if str(sample.peaks) != "nan":
                try:
                    sample = getFRiP(sample)
                except:
                    print("Record with FRiP value is empty or not found for sample %s" % sample.sampleName)
                    pass
        samples[sample.sampleName] = sample

    # write annotation sheet with statistics
    samples.T.to_csv(prj.sampleStats, index=False)

    print("Finished getting read statistics.")


def compare(args, prj):
    raise NotImplementedError


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
