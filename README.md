ChIP-seq pipelines
=========
ChIP-seq pipelines in Python
---------

[![Code Health](https://landscape.io/github/afrendeiro/chipseq-pipelines/master/landscape.svg?style=flat)](https://landscape.io/github/afrendeiro/chipseq-pipelines/master)

Tasks include raw sequencing quality control, adapter removal, read trimming, mapping, peak calling, peak annotation, target gene annotation, GO analysis and more. 



# Instalation
I recommend creating a virtual environment to install dependencies to running the pipeline.

First, make sure you work with python 2.7 or higher:
    
    module load python/2.7.6

Get virtualenv:

    curl -L -o virtualenv.py https://raw.githubusercontent.com/pypa/virtualenv/master/virtualenv.py

Create a virtual environment named `<name>`:

    python virtualenv.py <name>

Activate the environment:

    source <name>/bin/activate

Clone the repository

    git clone git@github.com:afrendeiro/chipseq-pipelines.git
    cd chipseq-pipelines

Install Python requirements:

    pip install -r requirements.txt

Good to go!

# Dependencies
It uses the following software, which is available in the CeMM cluster:
- Fastqc
- Piccard tools
- Bowtie2
- Samtools
- Bamtools

# Project structure
```
projectsroot
|__ projectname
    |__ runs
    |__ data
    |   |__ fastq
    |   |__ fastqc
    |   |__ raw
    |   |__ mapped
    |   |__ coverage
    |   |__ peaks
    |   |__ motifs
    |__ results
         |__ plots

htmlroot
|__ projectname
    |__ bigwig
```

JSON description [here](https://github.com/arendeiro/chipseq-pipelines/blob/master/projectPaths.json).

A thorough documentation of the folders will come in due time.
