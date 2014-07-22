# Concatenates all samples and make single file normalized file genome-wide
# Assumes coverage values from coverageBed are in column 4!!!

path <- "/sysdev/s3/share/data/oikopleura/chip-seq/mapped/"
outpath <- "/sysdev/s3/share/data/oikopleura/chip-seq/"
readsfile <- "/sysdev/s3/share/data/oikopleura/chip-seq/readcounts.tsv"
samples = paste(rep(c("TB", "D2", "D6-F", "D6-M"), each=4), c("IN", "E2F1", "E2F7", "POL2"), "R1", sep = "_")
# samples = paste(rep(c("TB", "D2", "D6-F", "D6-M"), each=8), rep(c("IN", "E2F1", "E2F7", "POL2"), each=2), rep(c("R1","R2"),2), sep = "_")
filenames = paste(path, samples, ".bwa.coverage.bed", sep = "")

raw = data.frame()

## For first sample:
### read in
sample1 <- read.table(filenames[1], sep="\t", header=FALSE)

### take genomic locations and coverage values for first sample
RAW = sample1[,c(1:4)]
colnames(RAW) <- c("chr", "start", "end", samples[1])

## add coverage values of other samples to table
for(i in 2:length(samples)) {
	print(paste("doing sample", samples[i], sep=" "))
	sampleTmp <- read.table(filenames[i], sep="\t", header=FALSE)
	RAW <- cbind(RAW, sampleTmp[4])
	colnames(RAW) <- c(colnames(RAW)[1:length(colnames(RAW))-1],samples[i])
}
# save as $ object
save(RAW, file=paste(outpath, "genome_coverage_raw.Rdata", sep=""))


# normalize the data

## get number of reads/sample
reads <- read.table(readsfile, sep="\t", header=FALSE)
colnames(reads) <- c("count","sample")

## find out which sample is "input"
inputs = grep("IN",colnames(RAW), ignore.case=TRUE)

NORM <- RAW[,1:3]

## select appropriate input for each sample
for (i in c(4:length(RAW))) {
	print(paste("doing sample", samples[i-3], sep=" "))
	if (i %in% inputs) {
		input <- i
		# get input name
		name <- samples[match(colnames(RAW[i]), samples)]
		# get total nº reads for that input
		totalIN <- reads$count[match(name, reads$sample)]
	} else {
		# get sample name
		name <- samples[match(colnames(RAW[i]), samples)]
		# normalize: log2, divide by input, normalize to total reads
		if (total / totalIN < 1) { 
			n <- log2((RAW[,i] + 0.001) / (RAW[,input] + 0.001) + (totalIN / total))
		} else {
			n <- log2((RAW[,i] + 0.001) / (RAW[,input] + 0.001) + (total / totalIN))
		}
		NORM <- cbind(NORM, n)
		colnames(NORM) <- c(colnames(NORM)[1:length(colnames(NORM))-1],colnames(RAW)[i])
	}
}

#save as R object
save(NORM, file=paste(outpath, "genome_coverage_norm.Rdata", sep=""))

# export as .bedgraph
for (i in 4:length(NORM)) {
	print(paste("doing sample", colnames(NORM[i]), sep=" "))
	write.table(NORM[c(1:3,i)], file=paste(outpath, "/normalized/", colnames(NORM[i]), ".bwa.norm.bedgraph", sep=""), sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)
}

# Plot correlations 
library(ggplot2)
library(reshape)

## Raw data
p = qplot(x=X1, y=X2, data=melt(cor(RAW[,c(4:length(RAW))])), fill=value, geom="tile") +
	xlab(NULL) + ylab(NULL) + theme_bw()
ggsave(filename = paste(outpath, "/plots/correlation.raw.png", sep=""), plot = p, height = 8, width = 15)
## Normalized
p = qplot(x=X1, y=X2, data=melt(cor(NORM[,c(4:length(NORM))])), fill=value, geom="tile") +
	xlab(NULL) + ylab(NULL) + theme_bw()
ggsave(filename = paste(outpath, "/plots/correlation.norm.png", sep=""), plot = p, height = 8, width = 15)
