Description of pipeline steps
=======

## Get static files required for analysis
get annotation, ChIP samples, etc...
get_files.sh

## Quality control

- QC (fastq)
fastqc_qc.sh

- Mapping (bwa) and read extension
bwa_mapping.sh

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
- Peak finding / ChIP enriched regions (peakzilla)
call_peaks.sh

- Peak genomic location
genome_distribution.sh
plot_genome_distribution.R

--------
- Binding motif (MEME)
find_motifs.sh

- Other TFs binding motifs

## TF binding context
- "Pseudo gene" enrichment
	TF ChIP signal enrichment in regions (TSS, TTS, all peaks - confirmation)

	Use "heatmap" model:
		- make "heatmap" model of 100 bp bins around region of interest (e.g. TSSs)
	TSS_windows.R

		- intersect normalized 200bp files with heatmap models on region of interest
	coverage_on_region.sh

	Make heatmap:
	plot_region_coverage.R

		- plot side by side (peaks, TSS, gene bodies?, TTS)

- Peaks enriched in Histone mods?
	- separated by genomic location (intergenic, TSS, intron)
	(interseption with ChIp enriched regions and statistical test vs random regions)

## TF target identification and characterization
- Target gene identification
	Decide on what parameters:
		- closeness (x Kb)
		- account for operons!
find_target_genes.sh

	Stage-specificity:
		Stage-specific targets set
		Common-stage targets set

- Target Validation:
	Expression correlation with TF binding
		TSS ChIP enrichment of E2F targets vs. non-targets/housekeeping genes
		Statistical test: expression levels targets vs. non-targets/housekeeping genes in appropriate stage
	plot_enrichment_targets_nontargets.R

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


