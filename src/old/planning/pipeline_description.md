Description of pipeline steps
=======

## Get static files required for analysis
get annotation, ChIP samples, etc...
get_annotation.sh

## Quality control

- QC (fastq)
fastqc_qc.sh

## Mapping and conversions
- Mapping (bwa) and read extension
bwa_mapping.sh

## Coverage, normalization and correlation
- Window coverage (200bp 100bp sliding windows)
	- Make windows in genome
	200bp_sliding_windows.sh

	- Calculate coverage on windows
	genome_coverage.sh

- Concatenation and input normalization
concatenate_and_normalize.R

- Sample correlation
sample_correlations.R

## TF binding characterization
- TF peak finding 
	
	This script calls peaks using peakzilla
	`call_peaks.sh`

- Peak statistics

	This script plots basic statistics on the peaks calculated by the peakfinder
	`peak_stats.R`

- Peak genomic location (genome distribution and distance to TSS (abs & rel))
	
	These scripts associate peaks and genomic locations, and plots comparisons with random locations
	`genome_distribution.sh`
	`plot_genome_distribution.R`

- Binding motif (MEME)

	This script calls peaks using MEME-chip
	`find_motifs.sh*`
	
	- Nº motifs per peak
	`meme2bed.py`
	# todo: `plot_motifs_peaks.R`

- Other TFs binding motifs


## TF binding context
- ChIP enrichment in regions (TSS, TTS, peaks)
	
	- This pipeline calculates a heatmap of ChIP enrichment around genomic regions of interest.

		1. Makes heatmap "model" of 100 bp bins around all regions of interest (e.g. TSSs)
		`TSS_windows.R`

		2. Computes ChIP coverage by intersecting ChIP files with heatmap models on region of interest
		`coverage_on_region.sh`

		3. Makes heatmap and line plot of averaged enrichment for all or a set of regions
		`plot_region_coverage.R`

- Peaks enriched in Histone mods?
	- separated by genomic location (intergenic, TSS, intron)
	(interception with ChIP enriched regions and statistical test vs random regions)

## TF target identification and characterization
- Target gene identification
	Decide on what parameters:
		- closeness (x Kb)
		- account for operons!
	`find_target_genes.sh`

	Stage-specificity:
		Stage-specific targets set
		Common-stage targets set

- Target Validation:
	Expression correlation with TF binding
		TSS ChIP enrichment of E2F targets vs. non-targets/housekeeping genes
		Statistical test: expression levels targets vs. non-targets/housekeeping genes in appropriate stage
	`plot_enrichment_targets_nontargets.R`
	`plot_expression_`


	Promoter activity correlation
		Heatmap Pol2, K27ac signals
		Statistical test: more enriched in targets vs. non-targets/housekeeping genes in appropriate stage

- Target's GO
target_genes_GO.R
	- All genes
	- Stage-specific
	- Common genes

- Specific target finding
	see suspected_targets.txt
	- cell cycle genes
	- histone genes

	Validation (qPCR)

## TF targets interspecies comparison
- Interspecies:
	- Orthologs
	- % conserved
	- Broad genomic position conserved
	- Specific position to TSS conserved (Distance)
	- What are the unconserved genes doing? (review GO and function)


