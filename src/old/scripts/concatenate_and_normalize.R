# Concatenates all samples and make single file normalized file genome-wide
# Assumes coverage values from coverageBed are in column 7!!!


## parse arguments (filenames)
path <- as.character(commandArgs(trailingOnly = TRUE))




require('gtools')
require('stringr')

path <- "/home/s3/afr/data/human/E2F7/genome_coverage/"
outpath <- "/home/s3/afr/data/human/E2F7/"
filenames = mixedsort(list.files(path))

## get a list of all samples to be processed
samples <- vector()
for (filename in filenames){
	sample = sub(pattern="[.].*",x=filename, replacement="")
	ifelse(test=sample %in% samples, yes="", no=samples <- append(samples,sample))
}

dataTable = data.frame()

## For first sample:
print(paste("doing sample", samples[1], sep=" "))
   
### read in
sample1 <- read.table(paste(path, filenames[1], sep=""), sep="\t", header=FALSE)

### take genomic locations and coverage values for first sample
dataTable = sample1[,c(1,2,3,7)]
colnames(dataTable) <- c("chr", "start", "end", samples[1])

## add coverage values of other samples to table
for(filename in filenames[-1]) {
	sampleName = sub(pattern="[.].*",x=filename, replacement="")
    print(paste("doing sample", sampleName, sep=" "))
	sampleTmp <- read.table(paste(path,filename, sep=""), sep="\t", header=FALSE)
	dataTable <- cbind(dataTable, sampleTmp[7])
	colnames(dataTable) <- c(colnames(dataTable)[1:length(colnames(dataTable))-1],sampleName)
}

save(dataTable, file=paste(outpath, "genome_coverage_raw.Rdata", sep=""))


# normalize the data

## find out which sample is "input"
input <- colnames(dataTable)[grep("input",colnames(dataTable), ignore.case=TRUE)]
inputIndex <- grep("input",colnames(dataTable), ignore.case=TRUE)

dataTable.norm <- dataTable[,1:3]

##select appropriate input for each sample!!!
for (sample in c(4:length(dataTable))) {
	if (sample != inputIndex) {
		norm <- log2(
					((dataTable[,sample]+1)/sum(dataTable[,sample]))
					/
					((dataTable[,input]+1)/sum(dataTable[,input]))
				)
		dataTable.norm <- cbind(dataTable.norm, norm)
        colnames(dataTable.norm) <- c(colnames(dataTable.norm)[1:length(colnames(dataTable.norm))-1],colnames(dataTable)[sample])
	}
}

#save as R object
save(dataTable.norm, file=paste(outpath, "genome_coverage_norm.Rdata", sep=""))
# export as .bed
write.table(dataTable.norm, file=paste(outpath, "genome_coverage_norm.bed", sep=""), sep="\t",quote=FALSE, col.names=TRUE)
