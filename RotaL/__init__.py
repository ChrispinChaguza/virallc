#!/usr/bin/env python

import os
import sys
from Bio import SeqIO
import glob
import argparse
import datetime
import urllib.request
import time
import random
import shutil
from RotaL.SeqLib import SeqLib
from RotaL.cmdOptions import cmdOptions

def main():
    cmdOptionsVals = cmdOptions.ReadCmdOptions()

    if (not shutil.which("blastn")) or (not shutil.which("blastn")):
        print("–Install NCBI BLAST software before using RotaL for rotavirus lineage classification")
    else:
        pass

    if (not shutil.which("mafft")):
        print("–Install MAFFT sequence alignment software before using RotaL for rotavirus lineage classification")
    else:
        pass

    dbpath = os.path.expanduser("~/db.rotavirus.lineages/")

    if (not os.path.exists(dbpath)) or (cmdOptionsVals['updatedb']):
        os.makedirs(dbpath,exist_ok=True)

        for urlFileName in SeqLib.dbLocation():
            print("Downloading and setting up the rotavirus lineage classification database...",end="\r")
            time.sleep(0.5)
            sys.stdout.flush()

            try:
                urlResponseObj = urllib.request.urlretrieve(urlFileName,dbpath+"/"+str(os.path.basename(urlFileName)))
            except:
                print("Error: Failed to setup the rotavirus lineage classification database... Check the internet connection")
                sys.exit()

    else:
        pass
            
    with open(cmdOptionsVals['outputLineageFile'],"w") as outputLineageHandle:
        if cmdOptionsVals['includeSequence']:
            lineageReportHeader = str(cmdOptionsVals['csvOutput']).join(["seqFastaFile", \
                                        "seqName", \
                                        "VirusName", \
                                        "segment", \
                                        "assignedLineageName", \
                                        "assignedLineageType", \
                                        "closestMatchedRefSeq", \
                                        "closestMatchIdentity", \
                                        "closestMatchCoverage", \
                                        "coverageFlag", \
                                        "assignedLineageInfo", \
                                        "nucSequence","\n"])
        else:
            lineageReportHeader = str(cmdOptionsVals['csvOutput']).join(["seqFastaFile", \
                                        "seqName", \
                                        "VirusName", \
                                        "segment", \
                                        "assignedLineageName", \
                                        "assignedLineageType", \
                                        "closestMatchedRefSeq", \
                                        "closestMatchIdentity", \
                                        "closestMatchCoverage", \
                                        "coverageFlag", \
                                        "assignedLineageInfo", \
                                        "\n"])

        outputLineageHandle.write(lineageReportHeader)

        if cmdOptionsVals['verboseOutput']:
            print(lineageReportHeader.strip())
        else:
            pass

        for seqFileName in cmdOptionsVals['inputSeqFiles']:
            if os.path.exists(seqFileName):
                for eachSeq in SeqIO.parse(seqFileName,"fasta"):
                    tmpSeqFile = "tmpR."+str(random.randrange(1,1000000))+"."+str(random.randrange(1,1000000))+"."+str(datetime.datetime.now()).replace(" ","")+".fas"
                    with open(tmpSeqFile,"w") as tmpSeqFileHandle:
                        tmpSeqFileHandle.write(">{seqName}\n{seq}\n".format(seqName=eachSeq.id,seq=str(eachSeq.seq)).replace("-",""))

                    if os.path.exists(tmpSeqFile):
                        assignedLineages = SeqLib.GetLineage(tmpSeqFile) 
                        os.remove(tmpSeqFile)

                        if not assignedLineages:
                            if cmdOptionsVals['verboseOutput']:
                                 print("Error: Sequence {seqID} not processed".format(seqID=eachSeq.id))
                            else:
                                pass

                            continue
                        else:
                            pass

                        for eachSeqLineage in assignedLineages.keys():
                            if cmdOptionsVals['includeSequence']:
                                lineageReport = assignedLineages[eachSeqLineage]
                                lineageReportSummary = str(cmdOptionsVals['csvOutput']).join(["{seqFastaFile}", \
                                        "{seqName}", \
                                        "{segment}", \
                                        "{VirusName}", \
                                        "{assignedLineageName}", \
                                        "{assignedLineageType}", \
                                        "{closestMatchedRefSeq}", \
                                        "{closestMatchIdentity}", \
                                        "{closestMatchCoverage}", \
                                        "{coverageFlag}", \
                                        "{assignedLineageInfo}", \
                                        "\n"]).format(seqFastaFile=seqFileName,
                                      seqName=lineageReport['seqName'],VirusName="Rotavirus A",
                                segment=lineageReport['segment'],
                                assignedLineageName=lineageReport['assignedLineageName'],
                                assignedLineageType=lineageReport['assignedLineageType'],
                                closestMatchedRefSeq=lineageReport['closestMatchedRefSeq'],
                                closestMatchIdentity=round(lineageReport['closestMatchIdentity'],2),
                                closestMatchCoverage=round(lineageReport['closestMatchCoverage'],2),
                                coverageFlag=lineageReport['coverageFlag'],
                                assignedLineageInfo=lineageReport['assignedLineageInfo'],
                                contactDetails=lineageReport['contactDetails'],
                                nucSequence=lineageReport['seq'],
                                sepChar=cmdOptionsVals['csvOutput'])
                            else:
                                lineageReport = assignedLineages[eachSeqLineage]
                                lineageReportSummary = str(cmdOptionsVals['csvOutput']).join(["{seqFastaFile}", \
                                        "{seqName}", \
                                        "{VirusName}", \
                                        "{segment}", \
                                        "{assignedLineageName}", \
                                        "{assignedLineageType}", \
                                        "{closestMatchedRefSeq}", \
                                        "{closestMatchIdentity}", \
                                        "{closestMatchCoverage}", \
                                        "{coverageFlag}", \
                                        "{assignedLineageInfo}", \
                                        "\n"]).format(seqFastaFile=seqFileName,
                                      seqName=lineageReport['seqName'],VirusName="Rotavirus A",
                                segment=lineageReport['segment'],
                                assignedLineageName=lineageReport['assignedLineageName'],
                                assignedLineageType=lineageReport['assignedLineageType'],
                                closestMatchedRefSeq=lineageReport['closestMatchedRefSeq'],
                                closestMatchIdentity=round(lineageReport['closestMatchIdentity'],2),
                                closestMatchCoverage=round(lineageReport['closestMatchCoverage'],2),
                                coverageFlag=lineageReport['coverageFlag'],
                                assignedLineageInfo=lineageReport['assignedLineageInfo'],
                                contactDetails=lineageReport['contactDetails'],
                                nucSequence=lineageReport['seq'],
                                sepChar=cmdOptionsVals['csvOutput'])

                            outputLineageHandle.write(lineageReportSummary)

                            if cmdOptionsVals['verboseOutput']:
                                print(lineageReportSummary.strip())
                            else:
                                pass
                    else:
                        pass

if __name__=="__main__":
    main()
