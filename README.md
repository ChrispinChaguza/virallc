# ViralLC: A package rapid assignment of viral lineage nomenclature

## Get the source code
```
git clone https://github.com/ChrispinChaguza/virallc.git
```

## About virallc

## Setup

The easist way to install ViralLC is using Conda (upcoming!).
```
conda install -c conda-forge virallc
conda install -c bioconda virallc
```

Another way to download ViralLC from GitHub and then manually setup the environment for the package 

```
git clone https://github.com/ChrispinChaguza/virallc.git
cd virallc
```

## Install the requirement packages

ViralLC requires several packages below which can be installed using Conda and Pip:

```
conda install -c bioconda mafft=7.526 -y
conda install -c conda-forge python=3.14.2 -y
conda install -c bioconda nextclade=3.18.1 -y
conda install -c conda-forge biopython=1.86 -y
conda install -c bioconda blast=2.16.0 -y
conda install -c conda-forge pandas=3.0.0 -y
conda install -c conda-forge networkx=3.6.1 -y
```

pip install gitdir==1.2.7

Alternatively, the packages can be installed as shown below.
```
#conda env export > environment.yml
conda env create -n virallc -f environment.yml 
```

Follow the instructions below to build and install ViralLC
```
python -m build 
pip install --force-reinstall dist/virallc-*.whl 
```

## Basic usage

### Setting up database

When running virallc for the first time, it will automatically download and setup the virallc database in the home directory (~/db.rotavirus.lineages/). However, if new lineages or sublineages have been assigned, you can redownload and update your local database as follows:
```
virallc database --setupdb

virallc database --updatedb

virallc database --version
```

### Assigning lineages

The simplest way to run virallc is to provide a single or separate multiple input FASTA file containing a single or multiple rotavirus A sequences.
```
virallc assign --sequences input.fasta --output report.tsv --database dbname
```
Or using the following shorthand options
```
virallc assign -s input.fasta -o report.tsv -d dbname
```

If you have multiple input FASTA files, you can run virallc as follows:
```
virallc assign --sequences input1.fasta input2.fasta input3.fasta --output report.tsv --database dbname
```

### Other options

In addition, the user can specify "--nucseq" or "-n" option to show nucleotide sequences of the query and closest matching reference sequences in the output report. To silence the output on the terminal, specify 
"--quiet" or "-q" option.

## Cite
To be updated!
