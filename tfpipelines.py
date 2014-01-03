#!/usr/bin/env python
"""
tfpipelines.py

Andre F. Rendeiro <afrendeiro@gmail.com>
2014

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 2 as
published by the Free Software Foundation.
"""

import sys, csv, os, logging
from optparse import OptionParser


def main():
    # option parser
    usage = 'python tfpipelines.py [OPTIONS] samples.tsv'
    parser = OptionParser(usage = usage)
    
    parser.add_option("-b", "--basedir",
    type="str", dest="basedir", default = os.path.abspath("."),
    help="workspace base directory")
    
    parser.add_option("-a", "--annotdir",
    type="str", dest="annotdir", default = os.path.abspath("./annotation"),
    help="directory with genome annotation")

    parser.add_option("-o", "--outdir",
    type="str", dest="outdir", default = os.path.abspath("./output"),
    help="output directory")

    parser.add_option("-l", "--log",
    type="str", dest="log", default = os.path.abspath("log.txt"),
    help="directory/filename to store log file to. Default (./log.txt)")

    # read arguments and options, assign to independent variables
    global options
    (options, args) = parser.parse_args()
    
    if len(args) > 2 or len(args) == 0:
        # return help message if argument number is incorrect
        print __doc__
        parser.print_help()
        sys.exit(0)

    # Log configuration
    logging.basicConfig(filename=options.log, level=logging.DEBUG, format='%(asctime)s %(message)s')
        # Here's the usage:
        #logging.debug('This message should go to the log file')
        #logging.info('So should this')
        #logging.warning('And this, too')

    # temporary list of suppported genomes
    global supportedSpecies
    supportedSpecies = ["human", "oikopleura"]

    # grab positional argument (samples file)
    samplesfile = args[0]

    # open samples file
    samples = readSamples(samplesfile)

    species = set()
    for sample in samples:
        species.add(sample[3])
    species = list(species)

    for specie in species:
        if (specie not in supportedSpecies):
            logging.debug('Tried to use unsupported "%s" specie' % specie)
            raise NotImplementedError("%s is not yet supported. Please see list of supported species in the documentation." %specie)


    # get annotation for all species in use
    for specie in species:
        getAnnotation(specie)
            

    # call mapping pipeline for each sample
    #os.system('sh %s/pipelines/asd.sh' % options.basedir)

    #write_some_files(computed, output.txt)

def readSamples(infile):
    """opens input file, read it and returns it"""
    try:
        with open(infile, 'r') as f:
            logging.info('Opened file "%s" in read mode' % infile)
            return [line.strip("\n").split("\t") for line in f]

    except IOError:
        logging.debug("Error opening '%s' file for reading samples" % infile)
        raise IOError("Sample file '%s' not openable or doesn't exist" % infile)


def getAnnotation(specie):
    """gets annotation"""
    logging.info('Trying to get annotation for "%s"' % specie)
    if (specie in supportedSpecies):
        os.system('./scripts/get_annotation.sh -s %s %s' % (specie, options.basedir))
    else:
        logging.debug('Tried to get annonation for unsupported "%s" specie' % specie)
        raise NotImplementedError("%s is not yet supported. Please see list of supported species in the documentation." %specie)


def writeFile(output, outfile):
    # write to file
    wr = csv.writer(outfile, delimiter='\t')
    wr.writerow(output)
    outfile.close()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("Program canceled by user!\n")
        sys.exit(0)