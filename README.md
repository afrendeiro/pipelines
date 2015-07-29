pipelines
=========

[![Build Status](https://travis-ci.org/afrendeiro/pipelines.svg?branch=master)](https://travis-ci.org/afrendeiro/pipelines) [![Code Health](https://landscape.io/github/afrendeiro/pipelines/master/landscape.svg?style=flat)](https://landscape.io/github/afrendeiro/pipelines/master)

## Instalation
Clone the repository

    git clone git@github.com:afrendeiro/pipelines.git
    cd pipelines

Install module:

    python setup.py install

or with `pip`:

    pip install .

If not root, add the `--user` to either command to install to your home directory.

Dependencies will be installed as well if these are not met.

You'll have now a `ngsProject` command in your `$PATH`. Type `ngsProject -h` to see all options.

Additionally, you'll also have a `.pipeline_config.yaml` file in your `$HOME` upon installation.

## Config file

The `.pipeline_config.yaml` in your `$HOME` contains configurations that are read by the pipeline objects.

A thorough description of the recognised configurations will soon follow.

## Sample annotation sheet
The sample annotation sheet is a csv file containing sample annotation. It is passed as argument to `ngsProject`.

Columns recognized are: `cellLine, numberCells, technique, ip, patient, patientID, sampleID, treatment, biologicalReplicate, technicalReplicate, experimentName, genome, unmappedBam`, but others might be present as well.

##### Required fields are: **`technique`**, **`genome`** and **`unmappedBam`**
