from setuptools import setup

setup(
    name="pipelines",
    version="0.1",
    description="Pipelines in Python.",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics"],
    keywords="bioinformatics, sequencing, ngs, ChIP-seq, ATAC-Seq",
    url="https://github.com/afrendeiro/chipseq-pipelines",
    author="Andre Rendeiro",
    author_email="arendeiro@cemm.oeaw.ac.at",
    license="GPL2",
    packages=["pipelines"],
    install_requires=["numpy", "pandas"],
    scripts=["pipelines/chipseq_pipeline"],
    include_package_data=True,
    zip_safe=False
)
