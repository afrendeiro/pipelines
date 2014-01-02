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
    usage = 'python script.py [OPTIONS] samples.txt'
    parser = OptionParser(usage=usage)
    
    parser.add_option("-b", "--basedir",
    type="str", dest="base", default='.',
    help="workspace base directory")
    
    parser.add_option("-a", "--annotdir",
    type="str", dest="annot", default='./annotation',
    help="directory with genome annotation")

    parser.add_option("-o", "--outdir",
    type="str", dest="out", default='./output',
    help="output directory")

    parser.add_option("-l", "--log",
    type="str", dest="log", default='log.txt',
    help="directory/filename to store log file to: default = log.txt")

    # read arguments and options
    (options, args) = parser.parse_args()

    if len(args) > 2 or len(args) == 0:
        # return help message if argument number is incorrect
        print __doc__
        parser.print_help()
        sys.exit(0)

    # Log configuration
    logging.basicConfig(filename=options.log,level=logging.DEBUG)
        # Here's the usage:
        #logging.debug('This message should go to the log file')
        #logging.info('So should this')
        #logging.warning('And this, too')

    # grab positional argument (samples file)
    samplesfile = args[0]

    # open samples file
    samples = readSamples(samplesfile)
    print samples
    os.system("ls")

    #write_some_files(computed, output.txt)

def readSamples(infile):
    # open input file, read it and returns it
    try:
        with open(infile, 'r') as f:
            return [line.strip("\n").split("\t") for line in f]

    except IOError:
        print("Sample file not openable or doesn't exist")

def writeFile(output, outfile):
    # write to file
    wr = csv.writer(outfile, delimiter='\t')
    wr.writerow(output)
    outfile.close()
    writeLog('Successfully wrote %d' % outfile, program.log)


def do_some_work(options, args):
    # does some computation
    for i, line in enumerate(infile):
        print "Doing...."
        row = line.split('\t')
        # do actual computation here and produce output
        output.append("something")
    writeLog('Done, computed %d' % computation, program.log)


class Log(object):
    """Log of operations"""
    def __init__(self, logfile):
        super(Log, self).__init__()
        logfile = logfile
        
    def writeLog(string, filename):
        # append line of string to log file of filename
        f = open(filename, 'a')
        f.write('%s %s\n' % (strftime("%H:%M:%S", localtime()), string))
        f.close()


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("Program canceled by user!\n")
        sys.exit(0)