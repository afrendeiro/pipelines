TF ChIP-seq pipelines
=========
Analysis of TF ChIP-seq data
---------
These are pipelines for analysis of ChIP-seq data of transcription factors. Perform all tasks from sample quality control to target gene annotation, GO analysis, enrichment correlations, etc, and will be the basis for my ChIP-seq project in Oikopleura.

One goal of this pipeline is to be reproducible and therefore runnable on any environment provided dependencies are satisfied.

This pipeline should work with TF ChIP-seq samples of any organism, and therefore initial specifications should be set.

Hardcoded paths, names, etc... are prohibited, everything should be passed as argument through the pipeline.

# Installation

# Dependencies

# Usage
tfpipelines.sh [-options] samplesfile

## Samples file
One sample per line with sample attributes tab-delimited: name, IP/input, filename, organism (see example_samples.txt)

# Documentation
see [Pipeline Description](planning/pipeline_description.md)