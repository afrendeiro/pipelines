for each region:
	for each bin:
		for each gene:
			average coverage through all 50bp windows
				concatenate in to sample/bin matrix


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

# find out the length of the files

## Row number 
dataRows <- length(unique(read.table(paste(path, filenames[1], sep=""), sep="\t", header=FALSE)[,5]))
dataCols <- length(mixedsort(filenames[grep(pattern=samples[1],x=filenames)]))

dataTable = as.data.frame(matrix(NA, ncol=dataCols+1, nrow=dataRows))
dataTable.mean = as.data.frame(matrix(NA, ncol=dataCols, nrow=length(samples)-1))

## For each sample:
for (sample in 1:length(samples)) {
	
	binFiles <- mixedsort(filenames[grep(pattern=samples[sample],x=filenames)])

	## Start with first bin
	### read in sample first bin file		
	bin1 <- read.table(paste(path, binFiles[1], sep=""), sep="\t", header=FALSE)
	bin1name <- gsub(x=str_extract(binFiles[1], "[.].*?[.]"), pattern="[.]",replacement="")

	### grab all bed entries for same gene (4)
	### (one could also filter genes (e.g. by expression))
	print(paste("doing sample '", samples[sample], "' bin ", bin1name, sep=""))

	tmpData <- data.frame(matrix(NA, ncol=2, nrow=dataRows))
	genes <- unique(bin1[,5])
	for (gene in 1:length(genes)){
		# find index of bins for a gene
		geneBins <- which(genes[gene] == bin1[,5])
		# average it
		geneBins.mean <- mean(bin1[c(geneBins), 4])
		# add to temp dataframe
		tmpData[gene,1:2] <- c(as.character(genes[gene]), mean(geneBins.mean))
	}

	### Add values for each gene in that bin to final table
	dataTable[,1:2] <- tmpData
	### Add mean values for all genes in that bin to final table
	dataTable.mean[1,1] <- mean(tmpData[,2])


	## Repeat process for all bins
	for(bin in 2:length(binFiles)) {
		### read in sample first bin file		
		binCur <- read.table(paste(path, binFiles[1], sep=""), sep="\t", header=FALSE)
		binCurname <- gsub(x=str_extract(binFiles[1], "[.].*?[.]"), pattern="[.]",replacement="")

		### grab all bed entries for same gene (4)
		### (one could also filter genes (e.g. by expression))
		print(paste("doing sample '", samples[sample], "' bin ", binCurname, sep=""))

		tmpData <- data.frame(matrix(NA, ncol=2, nrow=dataRows))
		genes <- unique(binCur[,5])
		for (gene in 1:length(genes)){
			# find index of bins for a gene
			geneBins <- which(genes[gene] == binCur[,5])
			# average it
			geneBins.mean <- mean(binCur[c(geneBins), 4])
			# add to temp dataframe
			tmpData[gene,1:2] <- c(genes[gene], mean(geneBins.mean))
		}

		### Add values for each gene in that bin to final table
		dataTable[,1:2] <- tmpData
		### Add mean values for all genes in that bin to final table
		dataTable.mean[1,1] <- mean(tmpData[,2])
	}
	
	# save file for each sample
	#write.table(dataTable, file=paste(outpath, samples[sample], ".TSS_coverage_norm.bed", sep=""), sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)
	save(dataTable, file=paste(outpath, samples[sample], ".TSS_coverage_norm.Rdata", sep=""))
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

png(paste(outpath, "TSS_line.png", sep=""))

    x<-unlist(colnames(dataTable.mean))
    y<-unlist(dataTable.mean[1,])

    plot(x, y, type="l")

dev.off()