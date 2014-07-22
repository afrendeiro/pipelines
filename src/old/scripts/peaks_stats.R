

path <- "/home/s3/afr/data/human/E2F7/TSS_coverage/"
outpath <- "~/output/human/E2F7/"


peaks = read.table("/home/s3/afr/output/human/E2F7/E2F7_Hela.peaks.tsv", sep="\t",header=TRUE,comment.char="/")


pdf(paste(outpath,"peaks_statistics.pdf",sep="/"))
	
	par(mfrow=c(2,3))

	names = colnames(peaks)

	for (col in 6:length(peaks)){
		boxplot(peaks[col], main = names[col])
	}

dev.off()