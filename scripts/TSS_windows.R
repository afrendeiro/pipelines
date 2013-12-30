
heatmap = function(binmid){
	#Function to 
	#Parses gene annotation and computes genomic windows around TSSs
	

	first <- binmid -100
	second <- binmid +100

	# gene annotation
	gns <- read.delim("/home/s3/afr/data/human/annotation/hg19/hg19.refseq.short.bed", sep="\t",header=F)
	colnames(gns) <- c("chr","start","end","ID","strand")

	# number of transcripts
	gns.singleTx <- as.character(gns$ID) %in% names(table(gns$ID)[table(gns$ID)==1])

	ID <- character(0)
	chr2 <- character(0)
	start <- integer(0)
	end <- integer(0)

	# do the plus strand genes
	gns.positive <- (gns[gns$strand=="+" & gns$end > 3101,])

	chrs <- unique(as.character(gns.positive$chr))



	for (chr in chrs){

		gns.chr <- gns.positive[gns.positive$chr==chr,]

		for (i in 1:nrow(gns.chr)){
			
		#	TSS2 <- c(TSS2,gns.chr$start[i])
		        TSS <- gns.chr$start[i]
		        ID <- c(ID,as.character(gns.chr$ID[i]))
		        chr2 <- c(chr2,chr)
		        start <- c(start,TSS + first)
		        end <- c(end,TSS + second)

		        
		}
	}


	# do the minus starnd genes
	gns.negative <- (gns[gns$strand=="-" & gns$start > 3101,])

	chrs <- unique(as.character(gns.negative$chr))
		
	for (chr in chrs){

		gns.chr <- gns.negative[gns.negative$chr==chr,]

		for (i in 1:nrow(gns.chr)){
			
			TSS <- gns.chr$end[i]
		    ID <- c(ID,as.character(gns.chr$ID[i]))
			chr2 <- c(chr2,chr)
		    start <- c(start,TSS - second)
		    end <- c(end,TSS - first)
		        
	    }
	}
			
	write.table(data.frame(chr2,start,end,ID),file=sprintf("heatmap_models_%d.txt",binmid), sep="\t",
		    row.names=FALSE,append=FALSE,quote=FALSE,col.names=FALSE) 	
	print(sprintf("heatmap_models_%d.txt",binmid))	            
}

binmids <- seq(-3000, 3000, 100)

for (bin in binmids){
	heatmap(bin)
}


q()
n

### Confirmations:
# Confirm all have same length
wc -l heatmap_models_TSS/heatmap_models_*

# Confirm there's no scientific notation
grep -i 'e' heatmap_models_*
	find ./ -type f -exec sed -i 's/e+07/0000000/g' {} \;

# Confirm start < end
for sample in -2000 -1900 -1800 -1700 -1600 -1500 -1400 -1300 -1200 -1100 -1000 -900 -800 -700 -600 -500 -400 -300 -200 -100 0 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000
do
awk 'int($2)< $3' heatmap_models_$sample.txt | wc -l
done

# Confirm start > 0
for sample in -2000 -1900 -1800 -1700 -1600 -1500 -1400 -1300 -1200 -1100 -1000 -900 -800 -700 -600 -500 -400 -300 -200 -100 0 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000
do
awk 'int($2)> 0' heatmap_models_$sample.txt | wc -l
done

