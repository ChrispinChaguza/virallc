#!/usr/bin/env python

import os
import sys
import argparse
from Bio import SeqIO
import multiprocessing
from itertools import combinations

def compareSeqs(i,j,k):
    seqMatch=0
    seqMisMatch=0

    for r,s in zip(str(i.seq).replace("-","N").replace("X","N").upper(),str(j.seq).replace("-","N").replace("X","N").upper()):
        if r==s and not (r=="N" or s=="N"):
            seqMatch=seqMatch+1
        else:
            if r=="N" or s=="N":
                continue
            else:
                seqMisMatch=seqMisMatch+1

    seqDist=round(seqMatch/(seqMatch+seqMisMatch)*100,4) if (seqMatch+seqMisMatch)!=0 else 0

    len1=round(len(i.seq),4)
    len2=round(len(j.seq),4)

    len1noN=round(len(str(i.seq).replace("-","N").replace("X","N").upper().replace("N","")),4)
    len2noN=round(len(str(j.seq).replace("-","N").replace("X","N").upper().replace("N","")),4)

    cov1=round(len1noN/len1*100,4) if len1!=0 else 0
    cov2=round(len2noN/len2*100,4) if len2!=0 else 0

    if False:
        if k:
            print(f"{i.id}\t{j.id}\t{len1}\t{len2}\t{len1noN}\t{len2noN}\t{cov1}\t{cov2}\t{seqMatch}\t{seqMisMatch}\t{seqDist}\n")
        else:
            pass

    return(f"{i.id}\t{j.id}\t{len1}\t{len2}\t{len1noN}\t{len2noN}\t{cov1}\t{cov2}\t{seqMatch}\t{seqMisMatch}\t{seqDist}")

def main():
    version = "1.0.1"
    options=argparse.ArgumentParser(sys.argv[0],
                usage=argparse.SUPPRESS,
                description='alnPairDist: A tool for calculating pairwise similarity of taxa in a multiple sequence alignment',
                prefix_chars='-',
                add_help=True,
                epilog='Written by Chrispin Chaguza, St Jude Children\'s Research Hospital, 2025')


    options.add_argument('--aln','-a',action='store',required=True,nargs=1,
                        metavar='alignment',dest='input',
                        help='Input multiple sequence alignment in fasta format')
    options.add_argument('--out','-o',action='store',required=False,nargs=1,
                        metavar='outfile',dest='outfile',default="aln.pairwise.similarity.tsv",
                        help='Output file containing pairwise similarity values')
    options.add_argument('--threads','-t',action='store',required=False,nargs=1,
                        metavar='threads',dest='threads',default=1,
                        help='Number of threads (default=1)')
    options.add_argument('--quiet','-q',action='store_false',default=True,
                        dest='verbose',help='Show progress')
    options.add_argument('--version','-v',action='store_true',default=False,
                        dest='version',help='Show software version')

    if len(sys.argv[:])>1: 
        if sys.argv[1]=="-v" or sys.argv[1]=="--version":
            print(f"alnPairDist {version}")
            sys.exit()
        else:
            options=options.parse_args(args=None if sys.argv[2:] else ['--help'])
    else:
        options=options.parse_args(args=None if sys.argv[2:] else ['--help'])

    cmdValues = {'inputSeqFiles': options.input[0:][0],
                 'outputFile': options.outfile[0:][0] if isinstance(options.outfile,list) else options.outfile,
                 'threads': int(options.threads[0:][0]) if isinstance(options.threads,list) else options.threads,
                 'verboseOutput': options.verbose,
                 'version': options.version}

    alignment={}

    for i in SeqIO.parse(cmdValues['inputSeqFiles'],"fasta"):
        alignment[i.id]=i

    threads=cmdValues['threads']

    fhandle=open(str(cmdValues['outputFile']),"w")

    if cmdValues['verboseOutput']:
        print("seq1\tseq2\tseq1.len\tseq2.len\tseq1.lenNoNs\tseq2.lenNoNs\tseq1.cov\tseq2.cov\tmatch\tmismatch\tpid\n")
    else:
        pass

    fhandle.write("seq1\tseq2\tseq1.len\tseq2.len\tseq1.lenNoNs\tseq2.lenNoNs\tseq1.cov\tseq2.cov\tmatch\tmismatch\tpid\n")

    tmpList = sorted(map(sorted, combinations(set([i for i in alignment.keys()]), 2)))

    with multiprocessing.Pool(processes=threads) as pool:
        args=[(alignment[m],alignment[n],cmdValues['verboseOutput']) for m,n in tmpList]
        #results=pool.starmap(compareSeqs, args)
        results=pool.starmap(compareSeqs, args)

    for result in results:
        if cmdValues['verboseOutput']:
            print(f"{result}")
        else:
            pass

        fhandle.write(f"{result}\n")

if __name__=="__main__":
    main()
