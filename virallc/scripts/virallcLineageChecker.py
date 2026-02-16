#!/usr/bin/env python

import os
import sys
import argparse
from pathlib import Path

def main():
    options=argparse.ArgumentParser(sys.argv[0],
                usage=argparse.SUPPRESS,
                description='virallcLineageChecker: A tool for checking novel virallc lineages',
                prefix_chars='-',
                add_help=True,
                epilog='Written by Chrispin Chaguza, St Jude Children\'s Research Hospital, 2025')

    options.add_argument('--lineage','-l',action='store',required=True,nargs=1,
                        metavar='lineage',dest='lineage',
                        help='Closest lineage inferred by virallc')
    options.add_argument('--name','-n',action='store',required=True,nargs=1,
                        metavar='seqname',dest='seqname',
                        help='Viral sequence or segment name in the database associated with lineage')
    options.add_argument('--db','-d',action='store',required=True,nargs=1,
                        metavar='database',dest='database',
                        help='virallc database to check')
    options.add_argument('--pid','-p',action='store',required=True,nargs=1,
                        metavar='nucidentity',dest='nucidentity',
                        help='Nucleotide identity inferred by virallc')
    options.add_argument('--cov','-c',action='store',required=True,nargs=1,
                        metavar='seqcoverage',dest='seqcoverage',
                        help='Sequence coverage inferred by virallc')
    options.add_argument('--dbfolder','-f',action='store',required=False,nargs=1,
                        metavar='databasedir',dest='databasedir',default='default',
                        help='Path to the virallc database folder')

    options=options.parse_args(args=None if sys.argv[2:] else ['--help'])

    cmdValues = {'lineage': options.lineage[0:][0],
                 'database': options.database[0:][0],
                 'seqname': options.seqname[0:][0],
                 'databasedir': options.databasedir[0:][0] if isinstance(options.databasedir,list) else options.databasedir,
                 'nucidentity': float(options.nucidentity[0:][0]),
                 'seqcoverage': float(options.seqcoverage[0:][0])}

    if cmdValues["databasedir"]=="default":
        cmdValues["databasedir"]=str(Path.home().joinpath("viraldb"))
        
        if os.path.exists(str(Path.home().joinpath("viraldb"))):
            pass
        else:
            print(f"Database directory {cmdValues['databasedir']} does not exist!")
            sys.exit()
    else:
        if os.path.exists(cmdValues["databasedir"]):
            pass
        else:
            print(f"Database directory {cmdValues['databasedir']} does not exist!")
            sys.exit()

    if os.path.exists(str(Path(cmdValues["databasedir"]).joinpath(cmdValues["database"]))):
        pass
    else:
        print(f"Database directory {str(Path(cmdValues['databasedir']).joinpath(cmdValues['database']))} does not exist!")
        sys.exit()

    if (cmdValues["seqcoverage"]<80 and cmdValues["seqcoverage"]>0 and cmdValues["seqcoverage"]<=100):
        print("Sequence coverage should be greater than 80% to assign novel lineages!")
        sys.exit()
    else:
        pass

    if (cmdValues["nucidentity"]<0 or cmdValues["nucidentity"]>100):
        print("Sequence identity should range 0 to 100%!")
        sys.exit()
    else:
        pass

    lineageClustersFile = str(Path(cmdValues['databasedir']).joinpath(cmdValues['database'],'lineages',f"{cmdValues['seqname']}.clusters.tsv"))
    if os.path.exists(lineageClustersFile):
        pass
    else:
        print(f"Lineage file {lineageClustersFile} does not exist!")
        sys.exit() 

    lineageThresholdFile = str(Path(cmdValues['databasedir']).joinpath(cmdValues['database'],'lineages','lineage.thresholds.txt'))
    if os.path.exists(lineageThresholdFile):
        pass
    else:
        print(f"Lineage threshold file {lineageThresholdFile} does not exist!")
        sys.exit()


    lineageList = str(cmdValues['lineage']).split('.')
    for i,j in enumerate(lineageList):
        if i>0:
            lineageList[i]=int(lineageList[i])
        else:
            pass

    lineagePrefix = lineageList[0]

    lineageThresholds = {}
    with open(lineageThresholdFile,"r") as fhandle:
        for i in fhandle:
            lineageThresholds[i.strip().split("\t")[0]]=[float(j) for j in i.strip().split("\t")[1:3]]

    lineageNames = {}
    tmpName = []

    if len(lineageList)==2:
        stem = lineageList[0]
    else:
        stem = ".".join([str(q) for q in lineageList[0:(len(lineageList)-1)]])

    if (cmdValues['nucidentity']>=lineageThresholds[cmdValues['seqname']][0]) and (cmdValues['nucidentity']<=lineageThresholds[cmdValues['seqname']][1]):
        with open(lineageClustersFile,"r") as fhandle1:
            for l in fhandle1:
                if len(lineageList)==2:
                    stem = ".".join([str(q) for q in lineageList[0:(len(lineageList))]])
                else:
                    stem = str(".".join([str(q) for q in lineageList[0:(len(lineageList)-1)]]))+"."
                if str(str(l).strip().split('\t')[1]).startswith(stem):
                    if len(lineageList)==2:
                        pass
                    else:
                        tmpName.append(int(str(str(l).strip().split('\t')[1]).replace(stem,'').split('.')[0]))
                else:
                    pass

        if str(stem).endswith('.'):
            print(f"Putative novel lineage designation: {stem}{sorted(tmpName)[-1]+1}")
        else:
            print(f"Putative novel lineage designation: {stem}.{sorted(tmpName)[-1]+1}")

    else:
        if (cmdValues['nucidentity']<lineageThresholds[cmdValues['seqname']][0]):
            stem = f"{lineageList[0]}."

            with open(lineageClustersFile,"r") as fhandle1:
                for l in fhandle1:
                    if str(str(l).strip().split('\t')[1]).startswith(stem):
                        tmpName.append(int(str(str(l).strip().split('\t')[1]).replace(stem,'').split('.')[0]))
                    else:
                        pass

            if str(stem).endswith('.'):
                print(f"Putative novel lineage designation: {stem}{sorted(tmpName)[-1]+1}")
            else:
                print(f"Putative novel lineage designation: {stem}.{sorted(tmpName)[-1]+1}")

        else:
            print(f"Putative novel lineage designation: {cmdValues['lineage']}")

if __name__=="__main__":
    main()



