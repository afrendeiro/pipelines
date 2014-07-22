
samples = list(
			c("output/oikopeura/ChIP-chip/E2F1_0.95.Female.location.bed",
			 "output/oikopeura/ChIP-chip/E2F1_0.95.Female.shuffled.location.bed"),
			c("output/oikopeura/ChIP-chip/E2F1_0.95.Male.location.bed",
			 "output/oikopeura/ChIP-chip/E2F1_0.95.Male.shuffled.location.bed"),
			c("output/oikopeura/ChIP-chip/E2F7_0.95.Female.location.bed",
			 "output/oikopeura/ChIP-chip/E2F7_0.95.Female.shuffled.location.bed"),
			c("output/oikopeura/ChIP-chip/E2F7_0.95.Male.location.bed",
			 "output/oikopeura/ChIP-chip/E2F7_0.95.Male.shuffled.location.bed")
			)

data <- data.frame()

samplePair<- samples[[1]]

for (samplePair in samples) {
	sample = samplePair[1]
	shuffled = samplePair[2]

	sampleName <- gsub(tail(strsplit(sample, split="/")[[1]], n=1), pattern=".[^.]*$", replacement="")
	sampleName <- gsub(sampleName, pattern=".location", replacement="")

	sampleData <- read.table(sample, sep="\t", header=FALSE)

	colnames(sampleData) <- c("chr", "start", "end", "score", "Location")
	sampleData$sample <- sampleName
	sampleData$type <- "Sample"

	sampleData.shuffled <- read.table(shuffled, sep="\t", header=FALSE)
	colnames(sampleData.shuffled) <- c("chr", "start", "end", "score", "Location")
	sampleData.shuffled$sample <- sampleName
	sampleData.shuffled$type <- "Shuffled"

	# append to data
	data <- rbind(data, sampleData, sampleData.shuffled)
}

data<-subset(data, !data$Location == ".")
save(data, file="output/oikopeura/ChIP-chip/E2F_genomicLocation.Rdata")

library("ggplot2")
library(extrafont)

p <- ggplot(data, aes(sample, fill=Location)) +
	geom_bar(position="fill") +
	facet_grid(. ~ type) +
	coord_flip() +
	xlab("Samples") +
	ylab("Count") +
	theme_bw() +
	theme(text=element_text(family="CMU Serif"), axis.text.x = element_text(angle = 50, hjust = 1))

ggsave(filename="output/oikopeura/ChIP-chip/E2F_genomicLocation.png",plot=p,height=7,width=10)

