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

        downloadBasePath = "http://127.0.0.1/"
        url = [downloadBasePath+"RVA.RefSeqs.db.mfa",
                downloadBasePath+"NSP1.db.final.fasta",
                downloadBasePath+"NSP2.db.final.fasta",
                downloadBasePath+"NSP3.db.final.fasta",
                downloadBasePath+"NSP4.db.final.fasta",
                downloadBasePath+"NSP5.db.final.fasta",
                downloadBasePath+"VP1.db.final.fasta",
                downloadBasePath+"VP2.db.final.fasta",
                downloadBasePath+"VP3.db.final.fasta",
                downloadBasePath+"VP4.db.final.fasta",
                downloadBasePath+"VP6.db.final.fasta",
                downloadBasePath+"VP7.db.final.fasta",
                downloadBasePath+"NSP1.db.final.aln",
                downloadBasePath+"NSP2.db.final.aln",
                downloadBasePath+"NSP3.db.final.aln",
                downloadBasePath+"NSP4.db.final.aln",
                downloadBasePath+"NSP5.db.final.aln",
                downloadBasePath+"VP1.db.final.aln",
                downloadBasePath+"VP2.db.final.aln",
                downloadBasePath+"VP3.db.final.aln",
                downloadBasePath+"VP4.db.final.aln",
                downloadBasePath+"VP6.db.final.aln",
                downloadBasePath+"VP7.db.final.aln",
                downloadBasePath+"NSP1.db.final.lineages.tsv",
                downloadBasePath+"NSP2.db.final.lineages.tsv",
                downloadBasePath+"NSP3.db.final.lineages.tsv",
                downloadBasePath+"NSP4.db.final.lineages.tsv",
                downloadBasePath+"NSP5.db.final.lineages.tsv",
                downloadBasePath+"VP1.db.final.lineages.tsv",
                downloadBasePath+"VP2.db.final.lineages.tsv",
                downloadBasePath+"VP3.db.final.lineages.tsv",
                downloadBasePath+"VP4.db.final.lineages.tsv",
                downloadBasePath+"VP6.db.final.lineages.tsv",
                downloadBasePath+"VP7.db.final.lineages.tsv"]

        for urlFileName in url:
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
            lineageReportHeader = "seqFastaFile{sepChar}" \
                                        "seqName{sepChar}" \
                                        "VirusName{sepChar}" \
                                        "segment{sepChar}" \
                                        "assignedLineageName{sepChar}" \
                                        "assignedLineageType{sepChar}" \
                                        "closestMatchedRefSeq{sepChar}" \
                                        "closestMatchIdentity{sepChar}" \
                                        "closestMatchCoverage{sepChar}" \
                                        "coverageFlag{sepChar}" \
                                        "assignedLineageInstr{sepChar}" \
                                        "nucSequence\n".format(sepChar=cmdOptionsVals['csvOutput'])
        else:
            lineageReportHeader = "seqFastaFile{sepChar}" \
                                        "seqName{sepChar}" \
                                        "VirusName{sepChar}" \
                                        "segment{sepChar}" \
                                        "assignedLineageName{sepChar}" \
                                        "assignedLineageType{sepChar}" \
                                        "closestMatchedRefSeq{sepChar}" \
                                        "closestMatchIdentity{sepChar}" \
                                        "closestMatchCoverage{sepChar}" \
                                        "coverageFlag{sepChar}" \
                                        "assignedLineageInstr{sepChar}" \
                                        "\n".format(sepChar=cmdOptionsVals['csvOutput'])

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
                                lineageReportSummary = "{seqFastaFile}{sepChar}" \
                                        "{seqName}{sepChar}" \
                                        "{segment}{sepChar}" \
                                        "{VirusName}{sepChar}" \
                                        "{assignedLineageName}{sepChar}" \
                                        "{assignedLineageType}{sepChar}" \
                                        "{closestMatchedRefSeq}{sepChar}" \
                                        "{closestMatchIdentity}{sepChar}" \
                                        "{closestMatchCoverage}{sepChar}" \
                                        "{coverageFlag}{sepChar}" \
                                        "{assignedLineageInstr}{sepChar}" \
                                        "\n".format(seqFastaFile=seqFileName,
                                      seqName=lineageReport['seqName'],VirusName="Rotavirus A",
                                segment=lineageReport['segment'],
                                assignedLineageName=lineageReport['assignedLineageName'],
                                assignedLineageType=lineageReport['assignedLineageType'],
                                closestMatchedRefSeq=lineageReport['closestMatchedRefSeq'],
                                closestMatchIdentity=round(lineageReport['closestMatchIdentity'],2),
                                closestMatchCoverage=round(lineageReport['closestMatchCoverage'],2),
                                coverageFlag=lineageReport['coverageFlag'],
                                assignedLineageInstr=lineageReport['assignedLineageInstr'],
                                contactDetails=lineageReport['contactDetails'],
                                nucSequence=lineageReport['seq'],
                                sepChar=cmdOptionsVals['csvOutput'])
                            else:
                                lineageReport = assignedLineages[eachSeqLineage]
                                lineageReportSummary = "{seqFastaFile}{sepChar}" \
                                        "{seqName}{sepChar}" \
                                        "{VirusName}{sepChar}" \
                                        "{segment}{sepChar}" \
                                        "{assignedLineageName}{sepChar}" \
                                        "{assignedLineageType}{sepChar}" \
                                        "{closestMatchedRefSeq}{sepChar}" \
                                        "{closestMatchIdentity}{sepChar}" \
                                        "{closestMatchCoverage}{sepChar}" \
                                        "{coverageFlag}{sepChar}" \
                                        "{assignedLineageInstr}{sepChar}" \
                                        "\n".format(seqFastaFile=seqFileName,
                                      seqName=lineageReport['seqName'],VirusName="Rotavirus A",
                                segment=lineageReport['segment'],
                                assignedLineageName=lineageReport['assignedLineageName'],
                                assignedLineageType=lineageReport['assignedLineageType'],
                                closestMatchedRefSeq=lineageReport['closestMatchedRefSeq'],
                                closestMatchIdentity=round(lineageReport['closestMatchIdentity'],2),
                                closestMatchCoverage=round(lineageReport['closestMatchCoverage'],2),
                                coverageFlag=lineageReport['coverageFlag'],
                                assignedLineageInstr=lineageReport['assignedLineageInstr'],
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
