#library("ChIPpeakAnno")
library(GenomicRanges)

path <- "/sysdev/s3/share/data/oikopleura/chip-seq/peaks/"
plotpath <- "/sysdev/s3/share/data/oikopleura/chip-seq/plots/"

samples = c(
			paste(rep(c("TB", "D2", "D6-F", "D6-M"), each=3), c("E2F1", "E2F7", "POL2"), "R1", sep = "_"),
			paste(paste(rep(c("TB", "D2", "D6-F", "D6-M"), each=3), c("E2F1", "E2F7", "POL2"), "R1", sep = "_"), "_shuffled", sep = "")
		)
# samples = paste(rep(c("TB", "D2", "D6-F", "D6-M"), each=6), rep(c("E2F1", "E2F7", "POL2"), each=2), rep(c("R1","R2"),2), sep = "_")
filenames = c(
			paste(path, paste(rep(c("TB", "D2", "D6-F", "D6-M"), each=3), c("E2F1", "E2F7", "POL2"), "R1", sep = "_"), ".bwa.peaks.annotation.tsv", sep = ""),
			paste(path, paste(rep(c("TB", "D2", "D6-F", "D6-M"), each=3), c("E2F1", "E2F7", "POL2"), "R1", sep = "_"), ".bwa.peaks.shuffled.annotation.tsv", sep = "")
		)

peaks = list()

# Make GRanges object for each sample
for(i in 1:length(samples)) {
	print(paste("Doing sample", samples[i], sep=" "))

	sample = read.table(filenames[i], sep="\t", header=FALSE)

	if (! grepl("D6", samples[i])) {
		colnames(sample) <- c("chr", "start", "end", "peakName", "peakSummit", "peakScore",
								"chipEnrich", "controlEnrich", "foldEnrich",
								"distScore", "FDR", "chrGL", "startGL", "endGL", "location",
								"overlapGL", "TSSchr", "TSSstart", "TSSend", "TSSstrand", "gene", "TSSdist"
							)
	} else {
		colnames(sample) <- c("chr", "start", "end", "peakName", "peakSummit", "peakScore",
								"chipEnrich", "controlEnrich", "foldEnrich",
								"distScore", "FDR", "chrGL", "startGL", "endGL", "location",
								"overlapGL", "TSSchr", "TSSstart", "TSSend", "TSSstrand", "gene", "TSSdist",
								"CHROMchr", "CHROMstart", "CHROMend", "CHROMstate", "CHROMdist"
							)
		sample$CHROMstate <- as.character(sample$CHROMstate)
	}
	# Add sample information
	sample$name <- samples[i]
	sample$stage <- strsplit(sample$name, "_")[[1]][1]
	sample$antibody <- strsplit(sample$name, "_")[[1]][2]
	sample$replicate <- strsplit(sample$name, "_")[[1]][3]
	sample$id <- paste(sample$name, sample$peakName, sep = "_")
	# Strings as character
	sample$location <- as.character(sample$location) ; sample$TSSstrand <- as.character(sample$TSSstrand); sample$gene <- as.character(sample$gene)

	# Compact in one line per peak again
	sample = aggregate(sample, by=list(sample$id), unique); sample$Group.1 <- NULL
	# Peaks annotated with more than one feature per type of annotation will have features comma-separated
	sample$location <- sapply(sample$location, paste, collapse = ",")
	sample$gene <- sapply(sample$gene, paste, collapse = ",")
	sample$TSSstart <- sapply(sample$TSSstart, paste, collapse = ",")
	sample$TSSstrand <- sapply(sample$TSSstrand, paste, collapse = ",")
	sample$TSSdist <- sapply(sample$TSSdist, paste, collapse = ",")

	if (! grepl("D6", samples[i])) {
		peaks [[ samples[i] ]] = with(sample,
			RangedData(
				IRanges(start = start, end = end, names = id),
				space = chr, stage = stage, antibody = antibody, replicate = replicate, summit = peakSummit, score = peakScore,
				chipEnrich = chipEnrich, controlEnrich = controlEnrich, foldEnrich = foldEnrich, FDR = FDR, location = location,
				gene = gene, TSSstart = TSSstart, TSSstrand = TSSstrand, TSSdist = TSSdist
			)
		)
	} else {
		sample$CHROMstate <- sapply(sample$CHROMstate, paste, collapse = ",")
		peaks [[ samples[i] ]] = with(sample,
			RangedData(
				IRanges(start = start, end = end, names = id),
				space = chr, stage = stage, antibody = antibody, replicate = replicate, summit = peakSummit, score = peakScore,
				chipEnrich = chipEnrich, controlEnrich = controlEnrich, foldEnrich = foldEnrich, FDR = FDR, location = location,
				gene = gene, TSSstart = TSSstart, TSSstrand = TSSstrand, TSSdist = TSSdist, CHROMstate = CHROMstate
			)
		)
	}
}

save(peaks, file = paste(path, "peak_annotation.Rdata", sep=""))


############ PLOTTING ###############
library(ggplot2)
library(reshape)
library(extrafont)


##### PEAKS ######

# Genomic location
df <- data.frame()

for (i in 1:length(peaks)) {
	if (grepl("shuffled", samples[i])) {type = "shuffled"} else {type = "data"}
	df <- rbind(df,
				data.frame(
					sample = samples[i],
					location = unlist(strsplit(peaks[[i]][["location"]],",")),
					type = type
				)
			)
}

p <- ggplot(df, aes(sample, fill = location)) +
	geom_bar(position = "fill") +
	facet_grid(. ~ type) +
	coord_flip() +
	xlab("Samples") +
	ylab("Frequency") +
	#guides(fill = guide_legend(ncol = 2)) +
	theme_bw() +
	theme(text = element_text(family = "CMU Serif"), axis.text.x = element_text(angle = 50, hjust = 1))

ggsave(filename = paste(plotpath, "/peak_genomic_locations.png", sep = ""), plot = p, height = 5, width = 7)


# Chromatin states

df <- data.frame()

for (i in 1:length(peaks)) {
	if (grepl("D6", samples[i])){
		if (grepl("shuffled", samples[i])) {type = "shuffled"} else {type = "data"}
		df <- rbind(df,
				data.frame(
					sample = samples[i],
					state = unlist(strsplit(peaks[[i]][["CHROMstate"]],",")),
					type = type
				)
			)
	}
}

p <- ggplot(df, aes(sample, fill = state)) +
	geom_bar(position = "fill") +
	coord_flip() +
	xlab("Samples") +
	ylab("Frequency") +
	guides(fill = guide_legend(ncol = 3)) +
	theme_bw() +
	theme(text = element_text(family = "CMU Serif"), axis.text.x = element_text(angle = 50, hjust = 1))

ggsave(filename = paste(plotpath, "/peak_chrom_states.png", sep = ""), plot = p, height = 5, width = 7)

# Distance to TSS
## Absolute

df <- data.frame()

for (i in 1:length(peaks)) {
	if (grepl("shuffled", samples[i])) {type = "shuffled"} else {type = "data"}
	df <- rbind(df,
			data.frame(
				sample = samples[i],
				start = as.numeric(start(peaks[[i]])),
				TSSstart = as.numeric(peaks[[i]][["TSSstart"]]),
				TSSstrand = peaks[[i]][["TSSstrand"]],
				TSSdist = as.numeric(peaks[[i]][["TSSdist"]]),
				type = type
			)
		)
}

# remove peaks without assigned TSS (".")
df <- df[df$TSSstart >= 0 ,]

p <- ggplot(df, aes(sample, TSSdist)) +
	geom_boxplot(outlier.shape = NA) +
	coord_flip(ylim = c(0,10000)) +
	#facet_grid(. ~ type) +
	ylab("Relative distance to TSS (bp)") +
	xlab("Samples") +
	theme_bw() +
	theme(text = element_text(family = "CMU Serif"), axis.text.x = element_text(angle = 50, hjust = 1))

ggsave(filename = paste(plotpath, "/peak_dist.abs_tss.png", sep = ""), plot = p, height = 5, width = 7)

## Relative
# Transform distance to negative if peak is upstream TSS
for (i in 1:nrow(df)) {
	if (df$TSSstrand == "+") {
		if (df$start[i] < df$TSSstart[i]) {
			df$TSSdist[i] <- -df$TSSdist[i]
		}
	} else {
		if (df$start[i] > df$TSSstart[i]) {
			df$TSSdist[i] <- -df$TSSdist[i]
		}
	}
}

p <- ggplot(df, aes(sample, TSSdist)) +
	geom_boxplot(outlier.shape = NA) +
	coord_flip(ylim = c(-5000,5000)) +
	#facet_grid(. ~ type) +
	ylab("Relative distance to TSS (bp)") +
	xlab("Samples") +
	theme_bw() +
	theme(text = element_text(family = "CMU Serif"), axis.text.x = element_text(angle = 50, hjust = 1))

ggsave(filename = paste(plotpath, "/peak_dist.rel_tss.png", sep = ""), plot = p, height = 5, width = 7)

##### GENES ######

## Raw data
p = qplot(x=X1, y=X2, data=melt(cor(RAW[,c(4:length(RAW))])), fill=value, geom="tile") +
	xlab(NULL) + ylab(NULL) + theme_bw()
ggsave(filename = paste(outpath, "/plots/correlation.raw.png", sep=""), plot = p, height = 8, width = 15)
## Normalized
p = qplot(x=X1, y=X2, data=melt(cor(NORM[,c(4:length(NORM))])), fill=value, geom="tile") +
	xlab(NULL) + ylab(NULL) + theme_bw()
ggsave(filename = paste(outpath, "/plots/correlation.norm.png", sep=""), plot = p, height = 8, width = 15)

