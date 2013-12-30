R

# Count
genLoc <- read.table("~/data/human/E2F7/E2F7_Hela.peaks.location.bed", sep="\t", header=FALSE)
colnames(genLoc) <- c("chr", "start", "end", "loc")

# Plot
genLoc <- read.table("~/data/human/E2F7/E2F7_Hela.peaks.location.bed", sep="\t", header=FALSE)
colnames(genLoc) <- c("chr", "start", "end", "loc")

peakLocs = data.frame(intergenic=0, TSS=0, utr5=0, exon=0, intron=0, utr3=0, TTS=0)

for (peak in genLoc$loc){
  if (peak == "."){
    peakLocs$intergenic <- peakLocs$intergenic+1
  } else {
    peak.loc = strsplit(x=as.character(peak),split="_")[[1]][3]
    if (peak.loc == "utr5"){
      peakLocs$utr5 <- peakLocs$utr5 + 1
    }
    if (peak.loc == "exon"){
      peakLocs$exon <- peakLocs$exon + 1
    }
    if (peak.loc == "intron"){
      peakLocs$intron <- peakLocs$intron + 1
    }
    if (peak.loc == "utr3"){
      peakLocs$utr3 <- peakLocs$utr3 + 1
    }
    if (peak.loc == "up"){
      peakLocs$TSS <- peakLocs$TSS + 1
    }
    if (peak.loc == "end"){
      peakLocs$TTS <- peakLocs$TTS + 1
    }
  }
}

peakLocs = t(peakLocs)

pie(peakLocs,labels=rownames(peakLocs))

# plot distance to TSS
TSS.distance <- read.table("~/data/human/E2F7/E2F7_Hela.peaks.TSS.bed", sep="\t",header=FALSE)
boxplot(TSS.distance$V5,outline=FALSE)

# plot relative 