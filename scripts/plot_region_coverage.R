# Plots heatmap and line plot of ChIP enrichment over given genomic locations
#
# Given a list of coverageBed files on a certain genomic regions, it:
## 1. Writes files of enrichment of each feature on that location
## 2. Writes files of mean feature enrichment for each sample
## 3. Plots heatmap of enrichment
## 4. Plots line
#
# Assumes coverage values from coverageBed are in column 8!!!
#
# André Rendeiro (afrendeiro@gmail.com) 

############# STUFF TO ADD ########################
# Arguments parsing

require('gtools')
require('stringr')

## parse arguments (filenames)
#path <- as.character(commandArgs(trailingOnly = TRUE))

path <- "/home/s3/afr/data/human/E2F7/TSS_coverage/"
outpath <- "~/output/human/E2F7/"
filenames = mixedsort(list.files(path))

## get a list of all samples to be processed
samples <- vector()
for (filename in filenames){
	sample = sub(pattern="[.].*",x=filename, replacement="")
	ifelse(test=sample %in% samples, yes="", no=samples <- append(samples,sample))
}


# find out which is input
input <- samples[grep("input",samples, ignore.case=TRUE)]
inputIndex <- grep("input",samples, ignore.case=TRUE)

# find out the length of the files
dataRows <- nrow(read.table(paste(path, filenames[1], sep=""), sep="\t", header=FALSE))
dataCols <- length(mixedsort(filenames[grep(pattern=samples[1],x=filenames)]))

dataTable = as.data.frame(matrix(NA, ncol=dataCols+1, nrow=dataRows))
dataTable.mean = as.data.frame(matrix(NA, ncol=dataCols, nrow=length(samples)-1))

## For each sample:
for (sample in 1:length(samples)) {
	if (sample != inputIndex) {
		
		binFiles <- mixedsort(filenames[grep(pattern=samples[sample],x=filenames)])
		INbinFiles <- mixedsort(filenames[grep(pattern=samples[inputIndex],x=filenames)])
		INreadCount <- as.numeric(read.table(paste(path, "../", samples[inputIndex], ".readcount.txt", sep=""), sep="\t", header=FALSE))

		### read in input first bin file
		INbin1 <- read.table(paste(path, INbinFiles[1], sep=""), sep="\t", header=FALSE)

		### read in sample first bin file		
		bin1 <- read.table(paste(path, binFiles[1], sep=""), sep="\t", header=FALSE)
		bin1name <- gsub(x=str_extract(binFiles[1], "[.].*?[.]"), pattern="[.]",replacement="")
		IPreadCount <- as.numeric(read.table(paste(path, "../", samples[sample], ".readcount.txt", sep=""), sep="\t", header=FALSE))

		### normalize
		### (one could also filter genes (e.g. by expression))
		print(paste("doing sample '", samples[sample], "' bin ", bin1name, sep=""))
		bin1Norm <- log2(
						((bin1[,8]+1)/IPreadCount)
						/
						((INbin1[,8]+1)/INreadCount)
					)
		dataTable[,1:2] = cbind(as.character(bin1[,4]),bin1Norm)
		# calculate mean for line plot ("pseudogene")
		dataTable.mean[sample,1] = mean(as.numeric(bin1Norm))

		colnames(dataTable)[1:2] <- c("gene", bin1name)
		colnames(dataTable.mean)[1] <- bin1name

		### add coverage values for remaining bins to table
		for(bin in 2:length(binFiles)) {
			# read sample bin
			binCur <- read.delim(paste(path, binFiles[bin], sep=""),  sep="\t", header=FALSE)
			binCurName <- gsub(x=str_extract(binFiles[bin], "[.].*?[.]"), pattern="[.]",replacement="")
			# read input bin
			INbinCur <- read.table(paste(path, INbinFiles[bin], sep=""), sep="\t", header=FALSE)
			# normalize
			print(paste("doing sample '", samples[sample], "' bin ", binCurName, sep=""))
			binCurNorm <- log2(
						((binCur[,8]+1)/IPreadCount)
						/
						((INbinCur[,8]+1)/INreadCount)
					)
			# Add data to dataTable and column names (bin identity)
			dataTable[,bin + 1] <- binCurNorm
			colnames(dataTable)[bin + 1] <- binCurName
			# calculate mean for line plot ("pseudogene")
			dataTable.mean[sample,bin] = mean(as.numeric(binCurNorm))
			colnames(dataTable.mean)[bin] <- binCurName
		}
		# save file for each sample
		#write.table(dataTable, file=paste(outpath, samples[sample], ".TSS_coverage_norm.bed", sep=""), sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)
		save(dataTable, file=paste(outpath, samples[sample], ".TSS_coverage_norm.Rdata", sep=""))
	}
}

# save file with all samples for line plot
#write.table(dataTable.mean, file=paste(outpath, "TSS_coverage_norm.mean.bed", sep=""), sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)
save(dataTable.mean, file=paste(outpath, "TSS_coverage_norm.mean.Rdata", sep=""))

# Heatmap 
load("~/output/human/E2F7/GSM810993_E2F7_Hela.TSS_coverage_norm.Rdata")
#dataTable<-read.table("~/output/human/E2F7/GSM810993_E2F7_Hela.TSS_coverage_norm.bed",sep="\t",header=TRUE)

pdf("heatmap.pdf")
heatmap(data.matrix(
	dataTable[,-1]),
	Colv=NA,
	keep.dendro=FALSE,
	scale="none"
	)
dev.off()

# Line Plot
#dataTable.mean<-read.table("~/output/human/E2F7/TSS_coverage_norm.mean.bed",sep="\t",header=TRUE)
load("~/output/human/E2F7/TSS_coverage_norm.mean.Rdata")
outpath <- "~/output/human/E2F7/"
pdf(paste(outpath, "TSS_line.pdf", sep=""))

    x<-unlist(colnames(dataTable.mean))
    y<-unlist(dataTable.mean[1,])

    plot(x, y, type="l")

dev.off()
