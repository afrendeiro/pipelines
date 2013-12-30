#!/usr/bin/env python
"""
<Shiny Program Name> <Program Version>

Andre F. Rendeiro <afrendeiro at gmail.com>
2012

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 2 as
published by the Free Software Foundation.
"""

import sys
import csv
import StringIO
import math
from optparse import OptionParser


def main():
    # option parser
    usage = 'python script.py [OPTIONS] file1.bed file2.bed > results.txt'
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--inputfile",
    type="str", dest="infile",
    help="file to open as input")

    # read arguments and options
    (options, args) = parser.parse_args()
    if len(args) > 1 or len(args) == 0:
        # return help message if argument number is incorrect
        print __doc__
        parser.print_help()
        sys.exit(0)
        infile = args[0]
    computed = do_some_work(options, args)
    #write_some_files(computed, output.txt)

def openFile(infile):
    # open input file, read it and returns it
    infile.open()
    return infile


def writeFile(output, outfile):
    # write to file
    wr = csv.writer(outfile, delimiter='\t')
    wr.writerow(output)
    outfile.close()

def do_some_work(options, args):
    # does some computation
    for i, line in enumerate(infile):
        print "Doing...."
        row = line.split('\t')
        # do actual computation here and produce output
        output.append("something")

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("Program canceled by user!\n")
        sys.exit(0)