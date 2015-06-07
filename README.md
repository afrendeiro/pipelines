ChIP-seq pipelines
=========
ChIP-seq pipelines in Python
---------

[![Code Health](https://landscape.io/github/afrendeiro/chipseq-pipelines/master/landscape.svg?style=flat)](https://landscape.io/github/afrendeiro/chipseq-pipelines/master)

Tasks include raw sequencing quality control, adapter removal, read trimming, mapping, peak calling, peak annotation, target gene annotation, GO analysis and more. 

# Instalation
I recommend creating a virtual environment to install dependencies to running the pipeline.
Get virtualenv:

    curl -L -o virtualenv.py https://raw.githubusercontent.com/pypa/virtualenv/master/virtualenv.py

Create a virtual environment named `<name>`:

    python virtualenv.py <name>

Activate the environment:

    source <name>/bin/activate

Clone the repository

    git clone git@github.com:afrendeiro/pipelines.git
    cd pipelines

Install module (it will install dependencies as well if these are not met):

    pip install .

You'll have now a `ngsProject` command in your `$PATH`. Type `ngsProject -h` to see all options.
