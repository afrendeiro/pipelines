# Plots line of ChIP enrichment over given genomic locations for specific sets of genes
#
# André Rendeiro (afrendeiro@gmail.com) 

############# STUFF TO ADD ########################
# Arguments parsing

path <- "/home/s3/afr/data/human/E2F7/TSS_coverage/"
outpath <- "~/output/human/E2F7/"

# Load normalized TSS enrichment data
load("~/output/human/E2F7/GSM810993_E2F7_Hela.TSS_coverage_norm.Rdata")

# Select target genes from universe
targets = unlist(read.table("/home/s3/afr/output/human/E2F7/E2F7_Hela.targets.txt",sep="\t",header=FALSE))
targetsRows = length(targets)

targetsTable = as.data.frame(matrix(NA, ncol=ncol(dataTable) - 1, nrow=targetsRows))
colnames(targetsTable) <- colnames(dataTable)[-1]

for (gene in 1:length(targets)){
	targetsTable[gene,] <- dataTable[dataTable$gene == targets[gene], -1]
}
targetsTable.mean = colMeans(data.matrix(targetsTable),na.rm=TRUE)

# Select control genes (housekeeping) from universe
controls = unlist(read.table("/home/s3/afr/data/human/annotation/HK_genes.txt",sep="\t",header=FALSE))
controlsRows = length(controls)

controlsTable = as.data.frame(matrix(NA, ncol=ncol(dataTable) - 1, nrow=controlsRows))
colnames(controlsTable) <- colnames(dataTable)[-1]

for (gene in 1:length(controls)){
	controlsTable[gene,] <- dataTable[dataTable$gene == controls[gene], -1]
}
controlsTable.mean = colMeans(data.matrix(controlsTable),na.rm=TRUE)

# Plot
outpath <- "~/output/human/E2F7/"

pdf(paste(outpath, "TSS_enrichment_targets_controls.pdf", sep=""))

    x<-unlist(colnames(targetsTable))
    y<-targetsTable.mean
    y2<-controlsTable.mean

	layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(4, 1))
		
		par(mai=rep(0.5, 4))

			plot_colors = c("red","gray")
			plot_names = c("target genes", "housekeeping genes")

		    plot(x, y, type="l", col="red",xlab="distance to TSS",ylab="fold enrichment")
	    	lines(x, y2, type="l", col="gray")

			boxplot(y, y2, col=plot_colors, names=plot_names)

			
		par(mai=c(0,0,0,0))
			plot.new()
			legend(x = "top",
					inset = 0,
	        		legend = plot_names,
	        		col = plot_colors,
	        		lwd=5,
	        		horiz = TRUE)
dev.off()

# Test if means are significantly different
wilcox.test(x=y, y=y2)

