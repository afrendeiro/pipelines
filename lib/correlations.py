#!/usr/env python

from argparse import ArgumentParser
import os
import re
from pybedtools import BedTool
import HTSeq
import numpy as np
import pandas as pd
import string
import itertools

import rpy2.robjects as robj  # for ggplot in R
import rpy2.robjects.pandas2ri  # for R dataframe conversion


def makeWindows(windowWidth, genome):
    """Generate 1kb windows genome-wide."""
    w = BedTool.window_maker(BedTool(), genome=genome, w=windowWidth)
    windows = dict()
    for interval in w:
        feature = HTSeq.GenomicInterval(
            interval.chrom,
            interval.start,
            interval.end
        )
        name = string.join(interval.fields, sep="_")
        windows[name] = feature

    return windows


def coverage(bam, intervals, fragmentsize, duplicates=False):
    """ Gets read coverage in bed regions, returns dict with region:count.
    bam - Bam object from HTSeq.BAM_Reader.
    intervals - dict with HTSeq.GenomicInterval objects as values.
    fragmentsize - integer.
    duplicates - boolean.
    """
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
    # Loop through TSSs, get coverage, append to dict
    cov = dict()

    for name, feature in intervals.iteritems():
        if feature.chrom not in chroms:
            continue
        count = 0

        # Fetch alignments in feature window
        for aln in bam[feature]:
            # check if duplicate
            if not duplicates and aln.pcr_or_optical_duplicate:
                continue

            # adjust fragment to size
            aln.iv.length = fragmentsize

            # add +1 to all positions overlapped by read within window
            count += 1

        # append feature profile to dict
        cov[name] = count
    return cov


def main(args):
    plotFunc = robj.r("""
        library(ggplot2)
        library(reshape2)

        function(df, path){
            # scatterplot
            pdf(path, height = 7, width = 7)

            par(pty = "s")

            smoothScatter(df[1], df[2], col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE, xaxt = 'n')
            text(par('usr')[1] + 1.8, par('usr')[4] - 0.5,
                bquote(R^2 == .(round(cor(A, B), 3))),
                cex = 1.6
            )
            title(xlab = colnames(df)[1],
                  ylab = colnames(df)[2],
                  outer = TRUE, cex.lab = 1.5)
            dev.off()
        }
    """)

    # Get sample names
    names = [re.sub(os.path.basename(bam), "\.bam", "") for bam in args.bamfiles]

    # Get genome-wide windows
    windows = makeWindows(args.window_width, args.genome)

    # Loop through all signals, compute coverage in bed regions, append to dict
    rawSignals = dict()
    for bam in xrange(len(args.bamfiles)):
        # Load bam
        bamfile = HTSeq.BAM_Reader(args.bamfiles[bam])
        # Get dataframe of bam coverage in bed regions, append to dict
        rawSignals[names[bam]] = coverage(bamfile, windows, args.fragment_size, args.duplicates)

    df = pd.DataFrame(rawSignals)

    # Normalize to library size
    dfNorm = df.apply(lambda x: np.log2(1 + (x / x.sum()) * 1000000))

    # pick samples pairwisely
    for sample1, sample2 in itertools.combinations(dfNorm.columns, 2):
        # pick two samples
        d = dfNorm[[sample1, sample2]]
        # convert the pandas dataframe to an R dataframe
        robj.pandas2ri.activate()
        df_R = robj.conversion.py2ri(d)
        # run the plot function on the dataframe
        plotFunc(df_R, os.path.join(args.plots_dir, sample1 + "_vs_" + sample2 + ".pdf"))


if __name__ == '__main__':
    # Parse command-line arguments
    parser = ArgumentParser(
        description='correlations.py',
        usage='python correlations.py [options] <directory> file1 [file2 ... fileN]'
    )

    # Global options
    # positional arguments
    parser.add_argument(dest='plots_dir', type=str, help='Directory to save plots to.')
    parser.add_argument('bamfiles', nargs='*', help='bamFiles')
    # optional arguments
    parser.add_argument('--duplicates', dest='duplicates', action='store_true')
    parser.add_argument('--window-width', dest='window_width', type=int, default=1000)
    parser.add_argument('--fragment-size', dest='fragment_size', type=int, default=50)
    parser.add_argument('--genome', dest='genome', type=str, default='hg19')

    args = parser.parse_args()

    main(args)
