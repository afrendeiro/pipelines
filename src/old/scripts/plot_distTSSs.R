# Since you already know some R I won't go much into details now

# Let's start by defining the location of the input file and the output folder
inputFile = "output/human/E2F7/GSM810993_E2F7_Hela.bwa.peaks.TSS.bed"
outPath = "output/human/E2F7/"
sample_name = "GSM810993_E2F7_Hela.bwa.peaks.TSS"

# Read in data and label columns
ChIPtoTSS <- read.table(inputFile, sep="\t",header=FALSE)
colnames(ChIPtoTSS) <- c("chr", "start", "end", ".",".","Gchr","Gstart", "Gend",".",".","strand", "distance")

# Remember we still have the ChIP-ER which didn't have a gene in the same scaffold in the file!
# They were annotated with distance "-1"
# We need to filter them:
ChIPtoTSS <- ChIPtoTSS[ChIPtoTSS$distance >= 0,]
# Here I'm simply asking for the rows ([row,]) which have distance >= 0
# in the column distance ($distance)

### PLOT
# plot distance to TSS
pdf(paste(outPath, sample_name,".peaks_absDist_TSS.pdf",sep=""))
	boxplot(ChIPtoTSS$distance,
		outline=FALSE,
		main="Absolute distance of peaks to TSSs",
		ylab="distance (bp)"
	)
dev.off()


# We can also calculate the distance relative to the gene
# (being 0 the TSS, upstream would be negative).

# Let's create the data.frame ChIPtoTSS.rel just to isolate variables.
ChIPtoTSS.rel <- ChIPtoTSS

# for each ChIP-ER/gene pair we ask if the gene is in the + ou - strand.
# if in the positive, we evaluate a further question:
# is the start of the ChIP-ER before (<) the start of the gene?
# if so, we replace the distance value with it's negative.
# The oposite is done for the negative strand:
for (row in 1:nrow(ChIPtoTSS.rel)) { # iterate through all ChIP-ER/gene pairs
	if (ChIPtoTSS.rel$strand[row] == "+") { # ask if gene is in + strand
		if (ChIPtoTSS.rel$start[row] < ChIPtoTSS.rel$Gstart[row]) { # if so, ask if ChIP-ER start < Gene start
			ChIPtoTSS.rel$distance[row] <- -ChIPtoTSS.rel$distance[row]	# if so, invert distance value
		}
	}
	if (ChIPtoTSS.rel$strand[row] == "-") { # ask if gene is in - strand
		if (ChIPtoTSS.rel$start[row] > ChIPtoTSS.rel$Gstart[row]) { # if so, ask if ChIP-ER start > Gene start
			ChIPtoTSS.rel$distance[row] <- -ChIPtoTSS.rel$distance[row]	# if so, invert distance value
		}
	}
}

# Plot the same way as above
pdf(paste(outPath, sample_name,".peaks_relDist_TSS.pdf",sep=""))
	boxplot(ChIPtoTSS.rel$distance,
		outline=FALSE,
		main="Relative distance of peaks to TSSs",
		ylab="distance (bp)"
	)
dev.off()


# Plot both in same
pdf(paste(outPath, sample_name,".peaks_Dist_TSS.pdf",sep=""))
	boxplot(
		ChIPtoTSS$distance,
		ChIPtoTSS.rel$distance,
		outline=FALSE,
		main="Distance of peaks to TSSs",
		ylab="distance (bp)",
		names=c("absolute","relative")
	)
dev.off()
