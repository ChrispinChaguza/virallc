#!/usr/bin/env python

import os
import sys
from Bio import SeqIO
from pathlib import Path

def main():

    for i in [Path(j) for j in sys.argv[1:]]:
        if i.name!=f"{str(i.stem)}.fasta":
            with open(f"{str(i.stem)}.fasta","w") as fhandle:
                    for j in SeqIO.parse(i,"fasta"):
                        sid=str(i).split(".")[0]
                        seq=j.seq

                        fhandle.write(f">{sid}\n{seq}\n")
        else:
            print(f"Input and output file have identical names (change the input file name)")
            pass

if __name__=="__main__":
    main()
