
infile_loc <- "output/oikopeura/ChIP-chip/ChIP-chip_Male_E2F7.peaks.location.bed"
infile_shuffled_loc <- "output/oikopeura/ChIP-chip/ChIP-chip_Male_E2F7.peaks.shuffled.location.bed"
infile_TSS <- "output/oikopeura/ChIP-chip/ChIP-chip_Male_E2F7.peaks.TSS.bed"
infile_shuffled_TSS <- "output/oikopeura/ChIP-chip/ChIP-chip_Male_E2F7.peaks.shuffled.TSS.bed"
outpath <- "output/oikopeura/ChIP-chip/"


##### GENOMIC LOCATIONS
sampleName <- gsub(tail(strsplit(infile_loc,split="/")[[1]], n=1),pattern=".[^.]*$", replacement="")

genLoc <- read.table(infile_loc, sep="\t", header=FALSE)
colnames(genLoc) <- c("chr", "start", "end", "score", "loc")

peakLocs = table(genLoc$loc)
peakLocs = t(peakLocs)
labels = paste(colnames(peakLocs), peakLocs, sep="\n")

#
genLocS <- read.table(infile_shuffled_loc, sep="\t", header=FALSE)
colnames(genLocS) <- c("chr", "start", "end", "score", "loc")

peakLocsS = table(genLocS$loc)
peakLocsS = t(peakLocsS)
labelsS = paste(colnames(peakLocsS), peakLocsS, sep="\n")

# plot peaks distribution
pdf(paste(outpath, "/", sampleName, ".pdf", sep = ""))

	par(mfrow=c(1,2))
	pie(peakLocs, labels = labels, main="peaks")
	pie(peakLocsS, labels = labelsS, main="shuffled")

dev.off()

# plot score per location
pdf(paste(outpath, "/", sampleName, ".score.pdf", sep = ""))
	boxplot(
		data = genLoc,
		score ~ loc,
		outline = FALSE,
		ylab = "peak score"
	)
dev.off()

##### TSSs
sampleName <- gsub(tail(strsplit(infile_TSS,split="/")[[1]], n=1),pattern=".[^.]*$", replacement="")

TSS.abs.distance <- read.table(infile_TSS, sep="\t",header=FALSE)
colnames(TSS.abs.distance) <- c("chr", "Pstart", "Gstart", "strand", "distance")

TSS.abs.distanceS <- read.table(infile_shuffled_TSS, sep="\t",header=FALSE)
colnames(TSS.abs.distanceS) <- c("chr", "Pstart", "Gstart", "strand", "distance")


# Test difference of means
wilcox = wilcox.test(x = TSS.abs.distance$distance, y = TSS.abs.distanceS$distance)

# plot distance to TSS
pdf(paste(outpath, "/", sampleName, ".abs_dist_TSS.pdf",sep=""))
	boxplot(
		TSS.abs.distance$distance,
		TSS.abs.distanceS$distance,
  		outline=FALSE,
  		main="Absolute distance of peaks to TSSs",
  		names=c("peaks", "shuffled"),
  		ylab="distance (bp)"
  	)

  	text(x=1.5, y=4000, labels = paste("p-value=", signif(wilcox$p.value, 2)))

dev.off()

TSS.distance <- TSS.abs.distance

for (row in 1:nrow(TSS.distance)) { 
	if (TSS.distance$strand[row] == "+") {
		if (TSS.distance$Pstart[row] < TSS.distance$Gstart[row]) {
			TSS.distance$distance[row] <- -TSS.distance$distance[row]
		}
	}
	if (TSS.distance$strand[row] == "-") {
		if (TSS.distance$Pstart[row] > TSS.distance$Gstart[row]) {
			TSS.distance$distance[row] <- -TSS.distance$distance[row]
		}
	}
}

TSS.distanceS <- TSS.abs.distanceS

for (row in 1:nrow(TSS.distanceS)) { 
	if (TSS.distanceS$strand[row] == "+") {
		if (TSS.distanceS$Pstart[row] < TSS.distanceS$Gstart[row]) {
			TSS.distanceS$distance[row] <- -TSS.distanceS$distance[row]
		}
	}
	if (TSS.distanceS$strand[row] == "-") {
		if (TSS.distanceS$Pstart[row] > TSS.distanceS$Gstart[row]) {
			TSS.distanceS$distance[row] <- -TSS.distanceS$distance[row]
		}
	}
}

# Test difference of means

wilcox = wilcox.test(x = TSS.distance$distance, y = TSS.distanceS$distance)

pdf(paste(outpath, "/", sampleName, ".rel_dist_TSS.pdf",sep=""))

	boxplot(
		TSS.distance$distance,
		TSS.distanceS$distance,
		outline=FALSE,
		main="Relative distance of peaks to TSSs",
  		names=c("peaks", "shuffled"),
		ylab="distance (bp)"
	)

	text(x=1.5, y=2000, labels=paste("p-value=", signif(wilcox$p.value, 2)))


dev.off()

