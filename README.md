# RotaL: A Python package for consistent assignment of rotavirus A lineages

## Get the source code
```
git clone https://github.com/ChrispinChaguza/rotavirus.lineage.classification.git
```

## About RotaL

## Basic usage

The simplest way to run RotaL is to provide a single or separate multiple input FASTA file containing a single or multiple rotavirus A sequences.
```
RotaL --sequences input.fasta --output report.tsv
```
Or using the following shorthand options
```
RotaL -s input.fasta -o report.tsv
```

If you have multiple input FASTA files, you can run RotaL as follows:
```
RotaL --sequences input1.fasta input2.fasta input3.fasta --output report.tsv
```

When running RotaL for the first time, it will automatically download and setup the rotavirus lineage classification database in the home directory (~/db.rotavirus.lineages/). However, if new lineages or sublineages have been assigned, you can redownload and update your local database as follows:
```
RotaL --sequences input.fasta --output report.tsv --updatedb
```

Or
```
RotaL -s input.fasta -o report.tsv -u
```

By default the output report containing the lineage calls will be provided in tab-separated value (TSV) file. To change the format of the output file to comma-separate value (CSV) format, specify "--csv" or "-c" options as follows:
```
RotaL --sequences input.fasta --output report.tsv --csv
```
Or
```
RotaL -s input.fasta -o report.tsv -c
```

In addition, the user can specify "--nucseq" or "-n" option to show nucleotide sequences of the query and closest matching reference sequences in the output report. To silence the output on the terminal, specify 
"--quiet" or "-q" option.

## Cite
Chrispin Chaguza, Chimwemwe Mhango, A. Duncan Steele, Carl D. Kirkwood, Jelle Matthijnssens, Francis E. Dennis, Martin M. Nyaga, Celeste M. Donato, and Khuzwayo C. Jere. ***RotaL: A Python package for consistent assignment of rotavirus A lineages***. Pending submission. 

Chrispin Chaguza, Celeste M. Donato, Chimwemwe Mhango, Arox W. Kamng’ona, Benjamin Kumwenda, Ernest Matambo, Daniel M. Njeru, A. Duncan Steele, Carl D. Kirkwood, Julius M. Mugweru, Jelle Matthijnssens, Rotavirus lineage designation consensus consortium, Francis E. Dennis, Martin M. Nyaga, and Khuzwayo C. Jere. ***Global rotavirus A evolution and genomic epidemiology revealed at high-resolution beyond genotypes with a new lineage classification nomenclature***. Pending submission.

