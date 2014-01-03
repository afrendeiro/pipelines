TF ChIP-seq pipelines
=========
Analysis of TF ChIP-seq data
---------
These are pipelines for analysis of ChIP-seq data of transcription factors. Perform all tasks from sample quality control to target gene annotation, GO analysis, enrichment correlations, etc, and will be the basis for my ChIP-seq project in Oikopleura.

One goal of this pipeline is to be reproducible and therefore runnable on any environment provided dependencies are satisfied.

This pipeline should work with TF ChIP-seq samples of multiple organisms, provided that initial specifications are set.

Hardcoded paths, names, etc... are prohibited, everything should be passed as argument through the pipeline.

# Installation

## With git:
	
	git clone git@github.com:afrendeiro/TF_ChIP-seq_pipelines.git
	cd TF_ChIP-seq_pipelines-master

## Manual download
	
	wget https://github.com/afrendeiro/TF_ChIP-seq_pipelines/archive/master.zip
	unzip master.zip
	cd TF_ChIP-seq_pipelines-master

# Dependencies
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