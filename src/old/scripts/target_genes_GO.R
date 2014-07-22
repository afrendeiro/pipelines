# Since you already know some R I won't go much into details now

# Let's start by defining the location of the input file and the output folder
inputFile = "output/oikopeura/ChIP-chip/ChIP-chip_Female_E2F7.location.bed"
outPath = "output/oikopeura/ChIP-chip/"
sample_name = "ChIP-chip_Female_E2F7"


CherGenes<-read.table("~/workspace/E2Fs/ChIP-chip/E2F7/Chipchip_chers_Male_E2F7_XC_rehybe_95_50.closestgenes.bed", sep="\t", header=F)

# distance to neigbour gene
colnames(CherGenes)[12]<-"distance"
mean(CherGenes$distance)
sd(CherGenes$distance)

pdf("distaneTSS.pdf")
	boxplot(CherGenes$distance,
		outline=FALSE,
  		main="Relative distance of peaks to TSSs",
  		ylab="distance (bp)"
  	)
dev.off()

# distance to gene that does not overlap with genes

# ChER overlapping more than one gene

# enriched GO terms
colnames(CherGenes)[9]<-"gene"

library("GSEABase")
library("GOstats")

# Background of all genes with GO annotations


if(file.exists("data//oikopleura//annotation//Oikopleura_GO_GeneSetCollection.Rdata")) {
  load("data//oikopleura//annotation//Oikopleura_GO_GeneSetCollection.Rdata")
} else {
  # create a GO frame annotation for Oikopleura
  
  library("GSEABase")
  library("GOstats")
  frame = read.table("data/oikopleura//annotation//Oiko_GO.go_id.txt",header=F)
  goframeData = data.frame(
    frame[,1],
    "ND",
    frame[,2])
  goFrame = GOFrame(goframeData, organism = "Oikopleura dioica")
  goAllFrame = GOAllFrame(goFrame)
  Oikopleura_GO_GeneSetCollection <- GeneSetCollection(goAllFrame, setType = GOCollection())
  save(Oikopleura_GO_GeneSetCollection, file="data/oikopleura//annotation/Oikopleura_GO_GeneSetCollection.Rdata")
}
universe = unique(unlist(geneIds(Oikopleura_GO_GeneSetCollection)))

# Set of genes to test
geneIDs = as.character(unique(unlist(CherGenes$gene)))
# don't forget to check if all peaks have genes (gene=".")
geneIDs = geneIDs[-1]

# Function to perform GO test:
GO.analysis<-function(genes, ontol, direction, universe){
  Test<-NA
  if(length(genes[genes %in% universe])>0){
    params <- GSEAGOHyperGParams(
      name = "GSEA based annot Params",
      geneSetCollection = Oikopleura_GO_GeneSetCollection,
      geneIds = genes,
      universeGeneIds = universe,
      ontology = ontol,
      pvalueCutoff = 0.01,
      conditional = FALSE,
      testDirection = direction)
    Test <- hyperGTest(params)
  }
  Test
}

## Tests for over-represented biological process, molecular function and cellular component terms:
#GO_BP<-GO.analysis(universe[1:1000], ontol="BP", direction="over", universe=as.list(universe))		# test using a subsample of universe
GO_BP<-GO.analysis(genes=geneIDs, ontol="BP", direction="over", universe=universe) #testing for overrepresentation here, but underrep also possible (or both)
GO_MF<-GO.analysis(genes=geneIDs, ontol="MF", direction="over", universe=universe)
GO_CC<-GO.analysis(genes=geneIDs, ontol="CC", direction="over", universe=universe)

## View results
summary(GO_BP)
summary(GO_MF)
summary(GO_CC)

# Use REViGO for visualisation
Output_for_REViGO<-function(GO_BP,GO_MF,GO_CC, cutoff){
  cutoff=0.01
  p.bp = pvalues(GO_BP)
  p.mf = pvalues(GO_MF)
  p.cc = pvalues(GO_CC)
  GO_data<-c(p.bp, p.mf, p.cc)
  GO_data_sig<-GO_data[GO_data<=cutoff]
  return(GO_data_sig)
}

GO_data_sig = Output_for_REViGO(GO_BP,GO_MF,GO_CC, cutoff=0.01)
write.table(GO_data_sig,file="GO_output_for_REViGO.txt",quote=FALSE,col.names=FALSE)

