TF ChIP-seq pipelines
=========
Compilation of pipelines and scripts for TF ChIP-seq data analysis
---------
These are pipelines for analysis of ChIP-seq data of transcription factors. Tasks include raw sample quality control, target gene annotation, GO analysis, enrichment correlations, and more. They are the basis for my ChIP-seq project in Oikopleura.

One goal of this repository is to increase reproducibility when reruning the same analysis.

These pipelines should work with TF ChIP-seq samples of multiple organisms, provided that anotation is provided.

# Installation

## With git:
	
	git clone git@github.com:afrendeiro/TF_ChIP-seq_pipelines.git
	cd TF_ChIP-seq_pipelines-master

# Dependencies
-
- Python
- R
	- Bioconductor

# Usage
	
	python tfpipelines.py [-options] samplesfile.tsv

Individual pipelines can be run individually, just add the '-h' option to see individual usage or read docs below.

IMPORTANT: some bash versions (non-GNU) will run slightly different than intended. The safest way to run the individual pipelines is by making them executable (`chmod +x <pipeline>`) and calling them directly.

### Samples configuration file
One sample per line with sample attributes tab-delimited: name, IP/input, filename, organism (see example samples.tsv)

# Documentation
see [Pipeline Description](planning/pipeline_description.md)