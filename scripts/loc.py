

infile = open("/home/afr/data/human/E2F7/E2F7_Hela.peaks.location.bed", "r")

exon=0
intron=0
utr5=0
utr3=0
promoter=0
TES=0

for i, line in enumerate(infile):
	peak = line.split('\t')[3]
	if "exon" in peak:
		exon += 1
	if "intron" in peak:
		intron += 1
	if "utr5" in peak:
		utr5 += 1
	if "utr3" in peak:
		utr3 += 1
	if "up" in peak:
		promoter += 1
	if "end" in peak:
		TES += 1

	

import numpy as np
import matplotlib.pyplot as plt

x = np.arange(0, 5, 0.1);
y = np.sin(x)
plt.plot(x, y)