# ViralLC: A package for rapid assignment of viral lineage nomenclature

## Getting the ViralLC source code
```
git clone https://github.com/ChrispinChaguza/virallc.git
```

## Setup ViralLC software on a local machine

### Installing ViralLC using Pip

The easist way to install the latest version of ViralLC is using Pip
```
pip install virallc
```

Here is a command to install a specific version of ViralLC using Pip (see the available versions on [PyPI](https://pypi.org/manage/project/virallc/releases/).
```
pip install virallc=={VERSION HERE}
```
After installing virallc using Pip, remember to install these dependencies (mafft, blast, and nextclade) manually using Conda!

### Installing ViralLC using Conda

Installation using Conda (upcoming!).
```
conda install -c conda-forge virallc
```
```
conda install -c bioconda virallc
```

### Installing ViralLC directly from Github

First, download ViralLC from GitHub and then manually setup the environment for the package 

```
git clone https://github.com/ChrispinChaguza/virallc.git
cd virallc
```

Second, manually install the required package dependencies (mafft, nextclade, biopython, blast, pandas, networkx, and gitdir) using Conda and Pip. If the installation fails, create a new Conda environment and then repeat the installation.

```
conda install -c bioconda mafft=7.526 -y
conda install -c conda-forge python=3.14.2 -y
conda install -c bioconda nextclade=3.18.1 -y
conda install -c conda-forge biopython=1.86 -y
conda install -c bioconda blast=2.16.0 -y
conda install -c conda-forge pandas=3.0.0 -y
conda install -c conda-forge networkx=3.6.1 -y
```
```
pip install gitdir==1.2.7
pip install build
```

Alternatively, the packages can be installed in a new Conda environment as shown below.
```
conda env create -n virallc -f environment.yml 
```

Follow the instructions below to build and install ViralLC
```
python -m build 
pip install --force-reinstall dist/{INSERT THE COMPILED SOFTWARE VERSION} 
```

## Basic usage

### Setting up database

When running virallc for the first time, it will automatically download and setup the virallc database in the home directory (~/viraldb/). However, if new lineages or sublineages have been assigned, you can redownload and update your local database as follows:

Here is a command to setup the database after installing the tool
```
virallc database --setupdb
```

```
virallc database -s
```

Here is a command to update a database (if corrupt, not found in the local machine, etc.)
```
virallc database --updatedb
```
```
virallc database -u
```

Here is a command to check version of the databases installed locally
```
virallc database --version
```
```
virallc database -v
```

### Assigning lineages

The simplest way to run virallc is to provide a single or separate multiple input FASTA file containing a single or multiple rotavirus A sequences.
```
virallc assign --in input.fasta --out report.tsv --db dbname
```
```
virallc assign -i input.fasta -o report.tsv -d dbname
```

To assign lineages to several sequences in a multi-FASTA file (each individual sequence represents a single strain):
```
virallc assign --in input1.fasta input2.fasta input3.fasta --out report.tsv --db dbname
```
```
virallc assign -i input1.fasta input2.fasta input3.fasta -o report.tsv -d dbname
```

To include the sequence in the output:
```
virallc assign --in input1.fasta input2.fasta input3.fasta --out report.tsv --db dbname --seq
```
```
virallc assign -i input1.fasta input2.fasta input3.fasta -o report.tsv -d dbname -s
```

To overwrite the output files:
```
virallc assign --in input1.fasta input2.fasta input3.fasta --out report.tsv --db dbname --seq --force
```
```
virallc assign -i input1.fasta input2.fasta input3.fasta -o report.tsv -d dbname -s -f
```

To assign lineages faster using more CPUs/threads:
```
virallc assign --in input1.fasta input2.fasta input3.fasta --out report.tsv --db dbname --seq --threads 10
```
```
virallc assign -i input1.fasta input2.fasta input3.fasta -o report.tsv -d dbname -s -t 10
```

To suppress the results on the terminal:
```
virallc assign --in input1.fasta input2.fasta input3.fasta --out report.tsv --db dbname --seq --threads 10 --quiet
```
```
virallc assign -i input1.fasta input2.fasta input3.fasta -o report.tsv -d dbname -s -t 10 -q
```

### Example dataset (Rotavirus A)

Here is a command to assign rotavirus A lineages to samples in *example* directory
```
virallc assign -i example.fna -o report.tsv -d RotavirusA
```

### Running ViralLC using Docker

Detailed information on installing ViralLC as a Docker container are available on [Docker Hub](https://hub.docker.com/r/chrispinchaguza/virallc).

### Full usage

General software options.
```
Contact: Chrispin Chaguza (Chrispin.Chaguza@STJUDE.ORG)

Usage:   virallc <command> [options]

Command: assign      assigns lineages to viral sequences
         database    setup, show, and update implemented databases
         version     prints program version
         citation    prints program citation information

Written by Chrispin Chaguza, St Jude Children's Research Hospital, 2025
```

Setup or update database options.

```
virallc: A tool for rapid assignment of virus lineages for given nomenclature

positional arguments:
  database

options:
  -h, --help      show this help message and exit
  --showdb, -p    Print list of implemented databases
  --setupdb, -s   Setup implemented databases
  --updatedb, -u  Update implemented databases
  --version, -v   Show database version

Written by Chrispin Chaguza, St Jude Children's Research Hospital, 2025
```

Assigning lineages to viral sequences.
```
virallc: A tool for rapid assignment of virus lineages for given nomenclature

positional arguments:
  assign

options:
  -h, --help            show this help message and exit
  --in, -i [query ...]  Input (multi-)fasta files to type (each contig is typed separately)
  --db, -d [refdb ...]  Specify viral database for the lineage classification (defaul=lineages.tsv)
  --out, -o outfile     Output file containing a summary of the assigned lineages
  --seq, -s             Show nucleotide sequence in the output
  --force, -f           Force overwrite output file
  --threads, -t threads
                        Number of threads (default=5)
  --quiet, -q           Show viral lineage assignment progress

Written by Chrispin Chaguza, St Jude Children's Research Hospital, 2025
```

Show software version.
```
virallc version
```

Show software citation information.
```
virallc citation
```

### Reference
```
Chrispin Chaguza et al., ViralLC: A package for rapid assignment of viral lineage nomenclature, GitHub, https://github.com/ChrispinChaguza/virallc.git
```

### Additional references

The following external software packages are run as part of the virallc software package.
```
gitdir, GitHub, https://github.com/sdushantha/gitdir
```
```
Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment search tool. J Mol Biol. 1990 Oct 5;215(3):403-10. doi: 10.1016/S0022-2836(05)80360-2. PMID: 2231712.
```
```
Katoh K, Misawa K, Kuma K, Miyata T. MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Res. 2002 Jul 15;30(14):3059-66. doi: 10.1093/nar/gkf436. PMID: 12136088; PMCID: PMC135756.
```
```
Aksamentov et al., (2021). Nextclade: clade assignment, mutation calling and quality control for viral genomes. Journal of Open Source Software, 6(67), 3773, https://doi.org/10.21105/joss.03773
```
