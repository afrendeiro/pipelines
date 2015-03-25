from setuptools import setup

setup(
    name="ChIPseq-pipelines",
    version="0.1",
    description="ChIP-seq pipelines in Python.",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics"],
    keywords="ATAC-Seq sequencing bioinformatics",
    url="https://github.com/afrendeiro/chipseq-pipelines",
    author="André Rendeiro",
    author_email="arendeiro@cemm.oeaw.ac.at",
    license="GPL2",
    packages=["chipseq_pipelines"],
    install_requires=["numpy", "pandas"],
    scripts=["bin/chipseq_pipelines"],
    include_package_data=True,
    zip_safe=False
)
