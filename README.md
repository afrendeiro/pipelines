ChIP-seq pipelines
=========
Compilation of pipelines and scripts for ChIP-seq data analysis
---------
Tasks include raw sequencing quality control, adapter removal, read trimming, mapping, peak calling, peak annotation, target gene annotation, GO analysis and more. 

# Dependencies
- bwa
- samtools (for now - to be replaced completely by sambamba)
- sambamba
- bedtools
- GNU parallel
- Python
- R
- R libraries (GenomicRanges, ggplot2, reshape, extrafont (optional))