#!/usr/bin/env python

"""
ChIP-seq pipeline

Workflow explained:
    - Project is created
    - Add Sample sheet to project (spawns next)
        - Samples are created and added to project (automatically)

In the process, stuff is checked:
    - project structure (created if not existing)
    - existance of csv sample sheet
    - existance of bam files from samples
    - read type of samples

:Example:

prj = Project("ngs")
prj.addSheet("sample_annotation.csv")
# that's it!

# explore!
prj.samples
prj.samples[0].mapped
prj.samples[3].nodupsshifted

prj.dirs.results
prj.sheet.to_csv(_os.path.join(prj.dirs.root, "sample_annotation.csv"))

# project options are read from the config file
# but can be changed on the fly:
prj = Project("ngs")
prj.config["mergetechnical"] = False
prj.addSheet("sample_annotation.csv")

"""

import os as _os
import pandas as _pd
import yaml as _yaml


class Paths(object):
    """
    A class to hold paths as attributes.
    """
    pass


class Project(object):
    """
    A class to model a Project.

    :param name: Project name.
    :type name: str

    Kwargs (will overule specified in config):
    :param parent: Path to where the project structure will be created.
    :type parent: str
    :param parenthtml: Path to where the project structure will be created.
    :type parenthtml: str

    :Example:

    prj = Project("ngs")
    prj = Project("ngs2", parent="/projects", parenthtml="/public_html")
    """
    def __init__(self, name, **kwargs):
        super(Project, self).__init__()
        # check it's a string
        self.name = name

        # Path structure
        self.dirs = Paths()

        # Read configuration file
        with open(_os.path.join(_os.path.expanduser("~"), "pipelines_config.yaml"), 'r') as handle:
            self.config = _yaml.load(handle)

        # If kwargs were passed, overule paths specified in the config with new ones.
        # parse kwargs
        for key in ('parent', 'parenthtml'):
            if key in kwargs:
                # check they're strings
                setattr(self.dirs, key, kwargs[key])
            else:
                setattr(self.dirs, key, self.config["paths"][key])

        # flow
        self.setProjectDirs()
        self.makeProjectDirs()
        self.setProjectPermissions()

        # samples
        self.samples = list()

    def __repr__(self):
        return "Project '%s'" % self.name

    def setProjectDirs(self):
        """
        Atributes directories for the project.
        """
        # check self.project_root exists and user has write access
        self.dirs.parent = _os.path.abspath(self.dirs.parent)
        if not _os.access(self.dirs.parent, _os.W_OK):
            raise IOError("%s does not exist, or user has no write access.\n\
            Use option '-r' '--project-root' to set a non-default project root path." % self.dirs.parent)

        # check self.html_root exists and user has write access
        self.dirs.parenthtml = _os.path.abspath(self.dirs.parenthtml)
        if not _os.access(self.dirs.parenthtml, _os.W_OK):
            raise IOError("%s does not exist, or user has no write access.\n\
            Use option '--html-root' to set a non-default html root path." % self.dirs.parenthtml)

        # Project directory structure
        self.dirs.root = _os.path.join(self.dirs.parent, self.name)
        # Runs
        self.dirs.runs = _os.path.join(self.dirs.root, "runs")
        self.dirs.pickles = _os.path.join(self.dirs.runs, "pickles")
        self.dirs.executables = _os.path.join(self.dirs.runs, "executables")
        self.dirs.logs = _os.path.join(self.dirs.runs, "logs")
        # Data
        self.dirs.data = _os.path.join(self.dirs.root, "data")
        self.dirs.fastq = _os.path.join(self.dirs.data, "fastq")
        self.dirs.fastqc = _os.path.join(self.dirs.data, "fastqc")
        self.dirs.raw = _os.path.join(self.dirs.data, "raw")
        self.dirs.mapped = _os.path.join(self.dirs.data, "mapped")
        self.dirs.coverage = _os.path.join(self.dirs.data, "coverage")
        self.dirs.peaks = _os.path.join(self.dirs.data, "peaks")
        self.dirs.motifs = _os.path.join(self.dirs.data, "motifs")
        # Results
        self.dirs.results = _os.path.join(self.dirs.root, "results")
        self.dirs.qc = _os.path.join(self.dirs.results, "qc")
        self.dirs.plots = _os.path.join(self.dirs.results, "plots")
        # Html structure
        self.dirs.html = _os.path.join(self.dirs.parenthtml, self.name)
        self.dirs.bigwig = _os.path.join(self.dirs.html, "bigwig")

        # Sample stats csv
        self.sampleStats = _os.path.join(self.dirs.root, self.name + ".sample_stats.csv")
        # Diffbind file
        self.diffBindCSV = _os.path.join(self.dirs.root, self.name + ".diffBind.csv")

    def makeProjectDirs(self):
        """
        Creates project directory structure if it doesn't exist.
        """
        for name, path in self.dirs.__dict__.items():
            if not _os.path.exists(path):
                _os.makedirs(path)

    def setProjectPermissions(self):
        """
        Makes the project's public_html folder executable.
        """
        for d in [self.dirs.parenthtml, self.dirs.html, self.dirs.bigwig]:
            try:
                _os.chmod(d, 0755)
            except OSError:
                # logger.error("cannot change folder's mode: %s" % d)
                continue

    def addSampleSheet(self, csv):
        """
        Build a `SampleSheet` object from a csv file and
        add it and its samples to the project.

        :param csv: Path to csv file.
        :type csv: str
        """
        # Make SampleSheet object
        self.sheet = SampleSheet(csv)

        # pair project and sheet
        self.sheet.project = self

        # Generate sample objects from annotation sheet
        self.sheet.getSamples()

        # Sample merging options:
        if self.config["options"]["mergetechnical"]:
            self.sheet.getBiologicalReplicates()
        if self.config["options"]["mergebiological"]:
            self.sheet.getMergedBiologicalReplicates()

        # Add samples to Project
        for sample in self.sheet.samples:
            # Check sample is from a suppoted genome
            if sample.genome not in self.config["genomes"]:
                raise TypeError("Sample's genome is not supported.")
            self.addSample(sample)
            sample.setFilePaths()

    def addSample(self, sample):
        """
        Adds a sample to the project's `samples`.
        """
        # Check sample is Sample object
        if not isinstance(sample, Sample):
            raise TypeError("Provided object is not a Sample object.")

        # Tie sample and project bilateraly
        sample.project = self
        # Append
        self.samples.append(sample)


class SampleSheet(object):
    """
    Class to model a sample annotation sheet.

    :param csv: Path to csv file.
    :type csv: str

    Kwargs (will overule specified in config):
    :param mergetechnical: Should technical replicates be merged to create biological replicate samples?
    :type mergetechnical: bool
    :param mergebiological: Should biological replicates be merged?
    :type mergebiological: bool

    :Example:

    prj = Project("ngs")
    sheet = SampleSheet("/projects/example/sheet.csv", prj)
    """
    def __init__(self, csv, **kwargs):

        super(SampleSheet, self).__init__()
        # TODO: checks on given args
        self.csv = csv

        # Read configuration file
        with open(_os.path.join(_os.path.expanduser("~"), "pipelines_config.yaml"), 'r') as handle:
            self.config = _yaml.load(handle)

        # Sample merging options:
        # If kwargs were passed, overule options specified in the config with new ones.
        # parse kwargs
        self.opts = dict()
        for key in ('mergetechnical', 'mergebiological'):
            if key in kwargs:
                # check they're strings
                self.opts[key] = kwargs[key]
            else:
                self.opts[key] = self.config["options"][key]

        self.samples = list()
        self._columns = ["cellLine", "numberCells", "technique", "ip", "patient",
                         "treatment", "biologicalReplicate", "technicalReplicate", "unmappedBam", "genome"]
        self.checkSheet()

    def __repr__(self):
        return "SampleSheet for project '%s'" % self.project

    def asDataFrame(self):
        """
        Returns a `pandas.DataFrame` representation of self.
        """
        return _pd.DataFrame(self.samples)

    def checkSheet(self):
        """
        Check if csv file exists and has all required columns.
        """
        try:
            self.df = _pd.read_csv(self.csv)
        except IOError("Given csv file couldn't be read.") as e:
            raise e

        missing = [col for col in self._columns if col not in self.df.columns]
        if len(missing) != 0:
            raise TypeError("Annotation sheet is missing columns: %s" % " ".join(missing))

    def getSamples(self):
        """
        Creates samples from annotation sheet and adds the samples to project.
        """
        for sample in range(len(self.df)):
            self.samples.append(Sample(self.df.ix[sample]))

    def getBiologicalReplicates(self):
        """
        Produces biological replicate samples from merged technical replicates.
        """
        # copy columns list
        attributes = self._columns[:]
        # ignore some fields in the annotation sheet
        attributes.pop(self._columns.index("unmappedBam"))
        attributes.pop(self._columns.index("technicalReplicate"))

        # get technical replicates -> biological replicates
        for key, values in self.df.groupby(attributes).groups.items():
            rep = self.df.ix[values][attributes].reset_index(drop=True).ix[0]
            if len(values) > 1:
                rep["technicalReplicate"] = 0
                rep["unmappedBam"] = self.df.ix[values]["unmappedBam"].tolist()
                # append biological replicate to samples
                self.samples.append(Sample(rep))

    def getMergedBiologicalReplicates(self):
        """
        Produces samples from merged biological replicates.
        """
        # copy columns list
        attributes = self._columns[:]
        # ignore some fields in the annotation sheet
        attributes.pop(self._columns.index("unmappedBam"))
        attributes.pop(self._columns.index("technicalReplicate"))
        attributes.pop(self._columns.index("biologicalReplicate"))

        # get biological replicates -> merged biological replicates
        for key, values in self.df.groupby(attributes).groups.items():
            rep = self.df.ix[values][attributes].reset_index(drop=True).ix[0]
            if len(values) > 1:
                # check samples in group are from different biological replicates
                if len(self.df.ix[self.df.groupby(attributes).groups[key]]['biologicalReplicate'].unique()) > 1:
                    rep["biologicalReplicate"] = 0
                    rep["technicalReplicate"] = 0
                    rep["unmappedBam"] = self.df.ix[values]["unmappedBam"].tolist()
                    # append merged biological replicate to samples
                    self.samples.append(Sample(rep))

    def to_csv(self, path):
        """
        Saves a csv annotation sheet from the samples.

        :param path: Path to csv file to be written.
        :type path: str

        :Example:

        sheet = SampleSheet("/projects/example/sheet.csv")
        sheet.to_csv("/projects/example/sheet2.csv")
        """
        df = _pd.DataFrame(self.samples)
        df = df[["sampleName"] + self._columns]
        # add extra columns with sample attributes
        df["mappedBam"] = [s.nodups for s in self.samples]
        df.to_csv(path, index=False)


class Sample(_pd.Series):
    """
    Class to model Samples basd on a pandas Series.

    :param series: Pandas `Series` object.
    :type series: pandas.Series

    :Example:

    prj = Project("ngs")
    sheet = SampleSheet("/projects/example/sheet.csv", prj)
    s1 = Sample(sheet.ix[0])
    """
    def __init__(self, series):
        from os.path import expanduser

        # Use _pd.Series object to have all sample attributes
        if not isinstance(series, _pd.Series):
            raise TypeError("Provided object is not a pandas Series.")
        super(Sample, self).__init__(series)

        self.checkNotEmpty()
        self.generateName()

        # Read configuration file
        with open(_os.path.join(expanduser("~"), "pipelines_config.yaml"), 'r') as handle:
            self.config = _yaml.load(handle)

        # check if sample is to be analysed with cuts
        self.tagmented = True if self.technique in self.config["tagment"] else False

        # Get track colour
        self.getTrackColour()

        # Get read type
        # self.getReadType()
        self.readType = "SE"
        self.paired = False

        # Get type of factor
        self.broad = True if self.technique in self.config["broadfactors"] else False
        self.histone = True if self.technique in self.config["histones"] else False

        # Only when sample is added to project, can paths be added.
        # The SampleSheet object, after beeing assigned to a project, will
        # call Sample.setFilePaths()

    def __repr__(self):
        return "Sample '%s'" % self.sampleName

    def checkNotEmpty(self):
        """
        Check if any important attribute is None.
        """
        req = ["genome", "biologicalReplicate", "technicalReplicate", "unmappedBam"]

        if not all([hasattr(self, attr) for attr in req]):
            raise ValueError("Required columns for sample do not exist.")

        if any([attr == "nan" for attr in req]):
            raise ValueError("Required values for sample are empty.")

    def generateName(self):
        """
        Generates a name for the sample by joining some of its attribute strings.
        """
        self.name = "_".join(
            str(x) for x in [
                self.cellLine, self.numberCells, self.technique,
                self.ip, self.patient, self.treatment,
                self.biologicalReplicate, self.technicalReplicate,
                self.genome
            ]
        )
        self.sampleName = self.name

    def asSeries(self):
        """
        Returns a `pandas.Series` object with all the sample's attributes.
        """
        return _pd.Series(self.__dict__)

    def getReadType(self, n=10):
        """
        Gets the read type (single, paired) and length of bam file.
        Returns tuple of (readType=string, readLength=int).

        :param n: Number of reads to read to determine read type. Default=10.
        :type n: int
        """
        import subprocess as sp
        from collections import Counter

        # for samples with multiple original bams, get only first
        if type(self.unmappedBam) == list:
            bam = self.unmappedBam[0]
        else:
            bam = self.unmappedBam

        # view reads
        p = sp.Popen(['samtools', 'view', bam], stdout=sp.PIPE)

        # Count paired alignments
        paired = 0
        readLength = Counter()
        while n > 0:
            line = p.stdout.next().split("\t")
            flag = int(line[1])
            readLength[len(line[9])] += 1
            if 1 & flag:  # check decimal flag contains 1 (paired)
                paired += 1
            n -= 1
        p.kill()

        # Get most abundant read length
        self.readLength = sorted(readLength)[-1]

        # If at least half is paired, consider paired end reads
        if paired > (n / 2):
            self.readType = "PE"
            self.paired = True
        else:
            self.readType = "SE"
            self.paired = False

    def setFilePaths(self):
        """
        Sets the paths of all files for this sample.
        """
        self.fastqc = self.project.dirs.fastqc

        fastq = _os.path.join(self.project.dirs.fastq, self.name)
        if self.readType == "SE":
            self.fastq = fastq + ".fastq"
        else:
            self.fastq1 = fastq + ".1.fastq"
            self.fastq2 = fastq + ".2.fastq"
            self.fastqUnpaired = fastq + ".unpaired.fastq"

        if self.readType == "SE":
            self.trimmed = fastq + ".trimmed.fastq"
        else:
            self.trimmed1 = fastq + ".1.trimmed.fastq"
            self.trimmed2 = fastq + ".2.trimmed.fastq"
            self.trimmed1Unpaired = fastq + ".1_unpaired.trimmed.fastq"
            self.trimmed2Unpaired = fastq + ".2_unpaired.trimmed.fastq"
        self.trimlog = _os.path.join(self.project.dirs.fastq, self.name + ".trimlog.txt")

        mapped = _os.path.join(self.project.dirs.mapped, self.name)
        self.mapped = mapped + ".trimmed.bowtie2.bam"
        self.alnRates = _os.path.join(self.project.dirs.results, self.name + ".alnRates.txt")
        self.alnMetrics = _os.path.join(self.project.dirs.results, self.name + ".alnMetrics.txt"),

        self.dups = mapped + ".trimmed.bowtie2.dups.bam"
        self.nodups = mapped + ".trimmed.bowtie2.nodups.bam"
        self.dupsMetrics = _os.path.join(self.project.dirs.results, self.name + ".duplicates.txt")

        # This will create additional bam files with reads shifted
        if self.tagmented:
            self.dupsshifted = mapped + ".trimmed.bowtie2.dups.shifted.bam"
            self.nodupsshifted = mapped + ".trimmed.bowtie2.nodups.shifted.bam"

        self.bigwig = _os.path.join(self.project.dirs.html, self.name + ".bigWig")
        self.trackURL = self.config["url"] + self.name + ".bigWig"
        self.coverage = _os.path.join(self.project.dirs.coverage, self.name + ".cov")

        # Analysis stuff
        self.peaksDir = _os.path.join(self.project.dirs.peaks, self.sampleName + "_peaks")
        self.peaks = _os.path.join(self.project.dirs.peaks, self.sampleName + "_peaks" + (".narrowPeak" if not self.broad else ".broadPeak"))
        self.frip = _os.path.join(self.project.dirs.results, self.sampleName + "_FRiP.txt")
        self.motifsDir = _os.path.join(self.project.dirs.motifs, self.sampleName)

        self.peaksMotifCentered = _os.path.join(self.peaksDir, self.sampleName + "_peaks.motifCentered.bed")
        self.peaksMotifAnnotated = _os.path.join(self.peaksDir, self.sampleName + "_peaks.motifAnnotated.bed")

    def getTrackColour(self):
        """
        Get a colour for a genome browser track based on the IP.
        """
        import random

        if self.ip in self.config["trackcolours"].keys():
            self.trackColour = self.config["trackcolours"][self.ip]
        else:
            if self.technique in ["ATAC", "ATACSEQ", "ATAC-SEQ"]:
                self.trackColour = self.config["trackcolours"]["ATAC"]
            elif self.technique in ["DNASE", "DNASESEQ", "DNASE-SEQ"]:
                self.trackColour = self.config["trackcolours"]["DNASE"]
            else:
                self.trackColour = random.sample(self.config["colourgradient"], 1)[0]  # pick one randomly
