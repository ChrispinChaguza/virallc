#!/usr/bin/env python

import os
import sys
import argparse
from Bio import SeqIO
import multiprocessing
from itertools import combinations

def compareSeqs(i,j):

    seqLen=round(len(i.seq),4)

    seqLennoN=round(len(str(i.seq).replace("-","N").replace("X","N").upper().replace("N","")),4)

    seqCov=round(seqLennoN/seqLen*100,4) if seqLen!=0 else 0

    return(f"{i.id}\t{seqLen}\t{seqLennoN}\t{seqCov}")

def main():
    options=argparse.ArgumentParser(sys.argv[0],
                usage=argparse.SUPPRESS,
                description='alnSeqCov: A tool for calculating sequence coverage of taxa in a multiple sequence alignment',
                prefix_chars='-',
                add_help=True,
                epilog='Written by Chrispin Chaguza, St Jude Children\'s Research Hospital, 2025')


    options.add_argument('--aln','-a',action='store',required=True,nargs=1,
                        metavar='alignment',dest='input',
                        help='Input multiple sequence alignment in fasta format')
    options.add_argument('--out','-o',action='store',required=False,nargs=1,
                        metavar='outfile',dest='outfile',default="aln.seq.cov.tsv",
                        help='Output file containing sequence coverage values')
    options.add_argument('--threads','-t',action='store',required=False,nargs=1,
                        metavar='threads',dest='threads',default=1,
                        help='Number of threads (default=1)')

    options=options.parse_args(args=None if sys.argv[2:] else ['--help'])


    cmdValues = {'inputSeqFiles': options.input[0:][0],
                 'outputFile': options.outfile[0:][0] if isinstance(options.outfile,list) else options.outfile,
                 'threads': int(options.threads[0:][0]) if isinstance(options.threads,list) else options.threads}


    alignment=[i for i in SeqIO.parse(cmdValues['inputSeqFiles'],"fasta")]

    threads=cmdValues['threads']

    fhandle=open(str(cmdValues['outputFile']),"w")

    fhandle.write("seq\tseq.len\tseq.lenNoNs\tseq.cov\n")

    with multiprocessing.Pool(processes=threads) as pool:
        args=[(m,m) for m in alignment]
        results=pool.starmap(compareSeqs, args)

    for result in results:
        fhandle.write(f"{result}\n")

if __name__=="__main__":
    main()
