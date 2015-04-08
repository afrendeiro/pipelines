ChIP-seq pipelines
=========
ChIP-seq pipelines in Python
---------

[![Code Health](https://landscape.io/github/afrendeiro/chipseq-pipelines/master/landscape.svg?style=flat)](https://landscape.io/github/afrendeiro/chipseq-pipelines/master)

Tasks include raw sequencing quality control, adapter removal, read trimming, mapping, peak calling, peak annotation, target gene annotation, GO analysis and more. 



# Installation
With pip it is as easy as:
    
    pip install git@github.com:afrendeiro/chipseq-pipelines.git

### Installing into a Python virtual environment
I recommend creating a virtual environment to install dependencies and to run.

Get virtualenv:
    
    sudo apt-get install python-virtualenv
or

    curl -L -o virtualenv.py https://raw.githubusercontent.com/pypa/virtualenv/master/virtualenv.py

Create a virtual environment named `<name>`:
    
    virtualenv <name>
or (in case you got it with curl as above)

    python virtualenv.py <name>

Activate the environment:

    source <name>/bin/activate

Install Python requirements in the virtual environment:

    pip install -r requirements.txt

Now install the package in the virtual environment:

    pip install git@github.com:afrendeiro/chipseq-pipelines.git

Good to go!

# Dependencies
It uses the following software:
- Fastqc
- Piccard tools
- Bowtie2
- Samtools
- Bamtools

make sure these executables are in your `$PATH`.


# Test
Run:
```sh
PROJECT_NAME=test
chipseq_pipeline --dry-run --no-checks preprocess $PROJECT_NAME examples/template_sample_annotation.csv
chipseq_pipeline --dry-run --no-checks analyse $PROJECT_NAME examples/template_sample_annotation.replicates.matchedControls.csv
chipseq_pipeline --dry-run --no-checks stats $PROJECT_NAME examples/template_sample_annotation.replicates.matchedControls.csv
chipseq_pipeline --dry-run --no-checks compare $PROJECT_NAME examples/template_sample_annotation.replicates.matchedControls.csv
```

for now you need to set these options as well until I start using a config parser:
`-r . --html-root . `

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
