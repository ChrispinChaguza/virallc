import sys
import os
from Bio import SeqIO
import random
import subprocess
import datetime
import random
from operator import itemgetter

dbpath = os.path.expanduser("~/db.rotavirus.lineages/")

def dbLocation():
    downloadBasePath = "https://raw.githubusercontent.com/ChrispinChaguza/RotaL/main/db.rotavirus.lineages/"
        
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

    return(url)

def LoadSeqAlignment(refSeqName):
    refAlignObj = {}

    for seqItem in SeqIO.parse(refSeqName,"fasta"):
        refAlignObj[seqItem.id] = str(seqItem.seq).upper()

    return(refAlignObj)

def NumSNPsFromSeqToOthers(seqAlign,seqName):
    seqObj = {}

    for seqItem in seqAlign.keys():
        if seqItem == seqName:
            continue

        seq1 = seqAlign[seqName].upper()
        seq2 = seqAlign[seqItem].upper()

        seqCoverage = 0
        seqMatches = 0

        for nucSeq1,nucSeq2 in zip(seq1,seq2):
            if (nucSeq1 not in ["n","N","-","X"]) and  (nucSeq2 not in ["n","N","-","X"]):
                seqCoverage = seqCoverage + 1

                if nucSeq1 == nucSeq2:
                    seqMatches = seqMatches + 1
                else:
                    continue
            else:
                continue

        seqObj[seqItem] = {"refSeqName": seqName, "seqCoverage": seqCoverage/len(seq1)*100,"seqMatches": seqMatches/seqCoverage*100}

    return(seqObj)

def NumSNPsBetweenSeqs(seqAlign,seqName1,seqName2):
    seqObj = {}
    seq1 = seqAlign[seqName1].upper()
    seq2 = seqAlign[seqName2].upper()
    seqCoverage = 0
    seqMatches = 0

    for nucSeq1,nucSeq2 in zip(seq1,seq2):
        if (nucSeq1 not in ["n","N","-","X"]) and  (nucSeq2 not in ["n","N","-","X"]):
            seqCoverage = seqCoverage + 1

            if nucSeq1 == nucSeq2:
                seqMatches = seqMatches + 1
            else:
                continue
        else:
            continue

    seqObj = {"seqName1": seqName1, "seqName2": seqName2, "seqCoverage": seqCoverage/len(seq1)*100,"seqMatches": seqMatches/seqCoverage*100}
    return(seqObj)

def GetClosestSeq(seqAlignName,refSeqName):
    seqAlign = LoadSeqAlignment(seqAlignName)
    snpDistObj = NumSNPsFromSeqToOthers(seqAlign,refSeqName)
    seqMinDist = {}

    for seqNum,seqName in enumerate(snpDistObj.keys()):
        if seqNum == 0:
            seqMinDist = {"closestSeqName": seqName, 'refSeqName': refSeqName, 
                    'seqMatches': snpDistObj[seqName]['seqMatches'], 'seqCoverage': snpDistObj[seqName]['seqCoverage']}
        else:
            if snpDistObj[seqName]['seqMatches'] > seqMinDist['seqMatches'] and seqName != refSeqName:
                seqMinDist = {"closestSeqName": seqName, 'refSeqName': refSeqName, 
                        'seqMatches': snpDistObj[seqName]['seqMatches'], 'seqCoverage': snpDistObj[seqName]['seqCoverage']}
            else:
                continue

    return(seqMinDist)

def AddSeqToAlignment(seqAlignFile,seqFile,threads=10):
    outAlignFile = "tmpR."+str(random.randrange(1,1000000))+"."+str(random.randrange(1,1000000))+"."+str(datetime.datetime.now()).replace(" ","")+".aln"

    mafftcmd = "mafft --quiet --auto --thread {threads} --keeplength --addfull {seqFile} \
            {seqAlignFile} > {outAlignFile}".format(seqFile=seqFile,seqAlignFile=seqAlignFile,outAlignFile=outAlignFile,threads=threads)
    subprocess.call(mafftcmd,shell=True,stdout=open(os.devnull,'w'),stderr=subprocess.STDOUT)

    if  os.path.exists(outAlignFile):
        tmpAlignSeqs = LoadSeqAlignment(outAlignFile)
        os.remove(outAlignFile)
        return(tmpAlignSeqs)
    else:
        return(False)

def WriteAlignment(seqObjDict,outputAlignFile):
    with open(outputAlignFile,"w") as alignFile:
        for seqName in seqObjDict.keys():
            alignFile.write(">{seqName}\n{seq}\n".format(seqName=seqName,seq=seqObjDict[seqName]))

def CompareSeqsBLAST(refSeqFile,seqFile,blastCoverageThresholdSegment):
    tmpBlastOutFile = "tmpR."+str(random.randrange(1,1000000))+"."+str(random.randrange(1,1000000))+"."+str(datetime.datetime.now()).replace(" ","")+".m8"
    blastcmd = "blastn -query {refSeqFile} -subject {seqFile} -outfmt \
            '6 qseqid sseqid qlen slen qstart qend length pident qcovs evalue qseq sseq' -out \
             {tmpBlastOutFile}".format(refSeqFile=refSeqFile,seqFile=seqFile,tmpBlastOutFile=tmpBlastOutFile)

    subprocess.call(blastcmd,shell=True,stdout=open(os.devnull,'w'),stderr=subprocess.STDOUT)

    blastMatch = {}

    if os.path.exists(tmpBlastOutFile):
        tmpFile = [str(record).strip() for record in open(tmpBlastOutFile,"r")]
    else:
        tmpFile = []

    if len(tmpFile) == 0:
        blastcmd = "tblastx -query {refSeqFile} -subject {seqFile} -outfmt \
                '6 qseqid sseqid qlen slen qstart qend length pident qcovs evalue qseq sseq' -out \
                 {tmpBlastOutFile}".format(refSeqFile=refSeqFile,seqFile=seqFile,tmpBlastOutFile=tmpBlastOutFile)
        subprocess.call(blastcmd,shell=True,stdout=open(os.devnull,'w'),stderr=subprocess.STDOUT)

        tmpFileRec = [str(record).strip() for record in open(tmpBlastOutFile,"r")]

        if len(tmpFileRec) == 0:
            return(False)
        else:
            pass
    else:
        pass

    tmpBlastMatchtmp = [str(record).strip().split("\t") for record in open(tmpBlastOutFile,"r")]

    for i,j in enumerate(tmpBlastMatchtmp):
        tmpBlastMatchtmp[i][6] = float(tmpBlastMatchtmp[i][6])

    tmpBlastMatch = sorted(tmpBlastMatchtmp,key=itemgetter(6), reverse=True)[0]

    os.remove(tmpBlastOutFile)

    blastMatch['qseqid'] = tmpBlastMatch[0]
    blastMatch['sseqid'] = tmpBlastMatch[1]
    blastMatch['qlen'] = tmpBlastMatch[2]
    blastMatch['slen'] = tmpBlastMatch[3]
    blastMatch['qstart'] = tmpBlastMatch[4]
    blastMatch['qend'] = tmpBlastMatch[5]
    blastMatch['length'] = tmpBlastMatch[6]
    blastMatch['pident'] = tmpBlastMatch[7]
    blastMatch['qcovs'] = tmpBlastMatch[8]
    blastMatch['evalue'] = tmpBlastMatch[9]
    blastMatch['qseq'] = tmpBlastMatch[10]
    blastMatch['sseq'] = tmpBlastMatch[11]

    if float(blastMatch['qcovs']) < blastCoverageThresholdSegment:
        return(False)
    else:
        return(blastMatch)

def GetAlignmentsForSegments(segmentName):
    rvaSegments = {"VP1": dbpath+"VP1.db.final.aln",
                    "VP2": dbpath+"VP2.db.final.aln",
                    "VP3": dbpath+"VP3.db.final.aln",
                    "VP4": dbpath+"VP4.db.final.aln",
                    "VP6": dbpath+"VP6.db.final.aln",
                    "VP7": dbpath+"VP7.db.final.aln",
                    "NSP1": dbpath+"NSP1.db.final.aln",
                    "NSP2": dbpath+"NSP2.db.final.aln",
                    "NSP3": dbpath+"NSP3.db.final.aln",
                    "NSP4": dbpath+"NSP4.db.final.aln",
                    "NSP5": dbpath+"NSP5.db.final.aln"}

    if segmentName in [i for i in rvaSegments.keys()]:
        return(rvaSegments[segmentName])
    else:
        return(False)

def GetSeqSeqForSegments(segmentName):
    rvaSegments = {"VP1": dbpath+"VP1.db.final.fasta",
                    "VP2": dbpath+"VP2.db.final.fasta",
                    "VP3": dbpath+"VP3.db.final.fasta",
                    "VP4": dbpath+"VP4.db.final.fasta",
                    "VP6": dbpath+"VP6.db.final.fasta",
                    "VP7": dbpath+"VP7.db.final.fasta",
                    "NSP1": dbpath+"NSP1.db.final.fasta",
                    "NSP2": dbpath+"NSP2.db.final.fasta",
                    "NSP3": dbpath+"NSP3.db.final.fasta",
                    "NSP4": dbpath+"NSP4.db.final.fasta",
                    "NSP5": dbpath+"NSP5.db.final.fasta"}

    if segmentName in [i for i in rvaSegments.keys()]:
        return(rvaSegments[segmentName])
    else:
        return(False)

def GetLineageThresholds(segmentName):
    thresholdPID = {"NSP1":[0.94,0.98],
                    "NSP2":[0.94,0.98],
                    "NSP3":[0.94,0.98],
                    "NSP4":[0.94,0.98],
                    "NSP5":[0.94,0.98],
                    "VP1":[0.93,0.97],
                    "VP2":[0.93,0.97],
                    "VP3":[0.93,0.97],
                    "VP4":[0.93,0.97],
                    "VP6":[0.93,0.97],
                    "VP7":[0.93,0.97]};

    return(thresholdPID[segmentName])

def GetRefLineageNames(segmentName):
    rvaSegments = {"VP1": dbpath+"VP1.db.final.lineages.tsv",
                    "VP2": dbpath+"VP2.db.final.lineages.tsv",
                    "VP3": dbpath+"VP3.db.final.lineages.tsv",
                    "VP4": dbpath+"VP4.db.final.lineages.tsv",
                    "VP6": dbpath+"VP6.db.final.lineages.tsv",
                    "VP7": dbpath+"VP7.db.final.lineages.tsv",
                    "NSP1": dbpath+"NSP1.db.final.lineages.tsv",
                    "NSP2": dbpath+"NSP2.db.final.lineages.tsv",
                    "NSP3": dbpath+"NSP3.db.final.lineages.tsv",
                    "NSP4": dbpath+"NSP4.db.final.lineages.tsv",
                    "NSP5": dbpath+"NSP5.db.final.lineages.tsv"}

    if segmentName in [i for i in rvaSegments.keys()]:
        return(rvaSegments[segmentName])
    else:
        return(False)

def GetRefSeqsSegmentsForBLAST():
    return(dbpath+"RVA.RefSeqs.db.mfa")

def assignLineageToSequence(segmentName,seqName):
    refLineageFile = GetRefLineageNames(segmentName)
    refLineageDB = {}

    for i in [str(i).strip().split("\t") for i in open(refLineageFile,"r")]:
        refLineageDB[i[0]] = i[1]
    if seqName in list(refLineageDB.keys()):
        return(refLineageDB[seqName])
    else:
        return(False)

def GetLineage(seqFastaFile):
    assignedLineageReport = {}
    blastCoverageThresholdSegment = 40

    for seqFastaItem in SeqIO.parse(seqFastaFile,"fasta"):
        tmpSeqName = seqFastaItem.id

        seqFile = "tmpR."+str(random.randrange(1,1000000))+"."+str(random.randrange(1,1000000))+"."+str(datetime.datetime.now()).replace(" ","")+".fas"

        with open(seqFile,"w") as tmpSeqHandle:
            tmpSeqHandle.write(">{seqID}____XX___\n{seq}\n".format(seqID=seqFastaItem.id,seq=seqFastaItem.seq))

        blastResults = CompareSeqsBLAST(GetRefSeqsSegmentsForBLAST(),seqFile,blastCoverageThresholdSegment)
        
        if not blastResults:
            return(False)
        else:
            pass

        if (not blastResults) or (float(blastResults['qcovs']) < blastCoverageThresholdSegment):
            assignedLineageReport = {"seqFastaFile": seqFastaFile, "seqName": tmpSeqName, "segment": 'NA', "assignedLineageName": 'NA',
                "assignedLineageType": 'NA', "closestMatchedRefSeq": 'NA',
                "closestMatchIdentity": 'NA', "closestMatchCoverage": 'NA',
                "coverageFlag": 'No match found', "assignedLineageInstr": 'NA', 
                "contactDetails": 'NA', "seq": str(seqFastaItem.seq)}
        else:
            refSegmentSeq = GetSeqSeqForSegments(blastResults['qseqid'])
            refSegmentALN = GetAlignmentsForSegments(blastResults['qseqid'])

            outSeqAlign = AddSeqToAlignment(refSegmentSeq,seqFile)
            del outSeqAlign[list(outSeqAlign.keys())[0]]

            #blastResults = CompareSeqsBLAST(GetRefSeqsSegmentsForBLAST(),seqFile,blastCoverageThresholdSegment)
            os.remove(seqFile)

            outSeqAlign.update(LoadSeqAlignment(refSegmentALN))
            outSeqAlignTMP = "tmpR."+str(datetime.datetime.now()).replace(" ","")+".aln"

            WriteAlignment(outSeqAlign,outSeqAlignTMP)

            closestMatchedRefSeq = GetClosestSeq(outSeqAlignTMP,blastResults['sseqid'])
            LineageName = assignLineageToSequence(blastResults['qseqid'],closestMatchedRefSeq['closestSeqName'])

            assignedLineageName = ""
            assignedLineageType = ""
            assignedLineageInfo = ""

            coverageFlag = ""
            contactDetails = "chrispin.chaguza@gmail.com"

            if closestMatchedRefSeq['seqMatches'] >= GetLineageThresholds(blastResults['qseqid'])[0]*100:
                if closestMatchedRefSeq['seqMatches'] >= GetLineageThresholds(blastResults['qseqid'])[1]*100:
                    assignedLineageName = LineageName 
                    assignedLineageType = "Known lineage"
                    assignedLineageInfo = "---"
                    assignedClosestLineage = LineageName
                else:
                    assignedLineageName = str(LineageName)+str(" [Novel?]")
                    assignedLineageType = "Potential new sublineage"
                    assignedClosestLineage = LineageName
                    assignedLineageInfo = "Submit full sequence to curators if coverage >85%"
            else:
                assignedLineageName = str(LineageName)+str(" [Novel?]")
                assignedLineageType = "Potential new lineage"
                assignedClosestLineage = LineageName
                assignedLineageInfo = "Submit full sequence to curators if coverage >85%"

            if closestMatchedRefSeq['seqCoverage'] >= 85:
                coverageFlag = "Excellent (lineage calls highly accurate)"
            elif closestMatchedRefSeq['seqCoverage'] >= 75:
                coverageFlag = "Good (lineage calls reasonably accurate)"
            elif closestMatchedRefSeq['seqCoverage'] >= 60:
                coverageFlag = "Fair (lineage calls maybe less accurate)"
            else:
                coverageFlag = "Poor (lineage calls likely inaccurate)"

            assignedLineageReport[tmpSeqName] = {"seqFastaFile": seqFastaFile, "seqName": tmpSeqName, "segment": blastResults['qseqid'], "assignedLineageName": assignedLineageName,
                "assignedLineageType": assignedLineageType, "closestMatchedRefSeq": closestMatchedRefSeq['closestSeqName'],
                "closestMatchIdentity": closestMatchedRefSeq['seqMatches'], "closestMatchCoverage": closestMatchedRefSeq['seqCoverage'],
                "coverageFlag": coverageFlag, "assignedLineageInfo": assignedLineageInfo, "contactDetails": contactDetails, "seq": str(seqFastaItem.seq).upper()}

            os.remove(outSeqAlignTMP)

    return(assignedLineageReport)














