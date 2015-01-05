#!/usr/bin/env python

"""
ChIP-seq pipeline 



"""

from argparse import ArgumentParser
import os
import logging
import pandas as pd

__author__ = "Andre Rendeiro"
__copyright__ = "Copyright 2014, Andre Rendeiro"
__credits__ = []
__license__ = "GPL3"
__version__ = "0.1"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"



def submit_job():
    pass

def bam2fastq():
    pass

def fastqc():
    pass

def trimm():
    pass

def bowtie2Map():
    pass

def shiftreads():
    pass

def shiftreads():
    pass

def markDuplicates():
    pass

def removeDuplicates():
    pass

def qc():
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
                                 "shifting", "markduplicates", "removeduplicates", "qc"],
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



    elif args.stage == "bam2fastq":
    elif args.stage == "fastqc":
    elif args.stage == "trimming":
    elif args.stage == "mapping":
    elif args.stage == "shifting":
    elif args.stage == "markduplicates":
    elif args.stage == "removeduplicates":
    elif args.stage == "qc":
