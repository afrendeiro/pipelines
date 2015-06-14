pipelines
=========

[![Code Health](https://landscape.io/github/afrendeiro/pipelines/master/landscape.svg?style=flat)](https://landscape.io/github/afrendeiro/pipelines/master)

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

Install module:

    python setup.py install

with `pip`, dependencies will be installed as well if these are not met:

    pip install .

You'll have now a `ngsProject` command in your `$PATH`. Type `ngsProject -h` to see all options.
