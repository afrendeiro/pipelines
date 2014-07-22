#!/bin/bash

bedtools makewindows -g ~/data/human/annotation/hg19/hg19_ChrSizes.tsv -w 200 -s 100 > ~/data/human/annotation/hg19/200bp_100bpsliding_windows.bed