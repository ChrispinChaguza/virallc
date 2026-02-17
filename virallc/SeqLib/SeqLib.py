import sys
import os
from Bio import SeqIO
import random
import subprocess
import datetime
from operator import itemgetter
import argparse
import shutil
import urllib.request
import concurrent.futures
import fsspec
from pathlib import Path

dbpath = str(Path.home().joinpath("viraldb"))
tmpdirpath = f"tmpR.virallc.{random.randrange(1,1000000)}.{datetime.datetime.now()}".\
                 replace(':','').replace(' ','.').replace(',','.')
refdbName = ""
version = "1.0.14"
name = "Chrispin Chaguza"
email = "Chrispin.Chaguza@STJUDE.ORG"

def ReadCmdOptions():
    if len(sys.argv)==1 or not (sys.argv[1] in ["assign","database","version","citation"]):
        print(f"Program: virallc (program for rapid viral lineage assignment)\n"\
               "Version: {version}\n"\
               "Contact: {name} ({email})\n\n"\
               "Usage:   virallc <command> [options]\n\n"\
               "Command: assign      assigns lineages to viral sequences\n"\
               "         database    setup, show, and update implemented databases\n"\
               "         version     prints program version\n"\
               "         citation    prints program citation information\n\n"\
               "Written by Chrispin Chaguza, St Jude Children\'s Research Hospital, 2025")

        if len(sys.argv)>1 and (not (sys.argv[1] in ["assign","database","version","citation"])):
            print(f"Unrecognised program option \"{sys.argv[1]}\"")
            sys.exit()
        else:
            sys.exit()

    if sys.argv[1]=="assign":
        options=argparse.ArgumentParser(sys.argv[0],
                    usage=argparse.SUPPRESS,
                    description='virallc: A tool for rapid assignment of virus lineages for given nomenclature',
                    prefix_chars='-',
                    add_help=True,
                    epilog='Written by Chrispin Chaguza, St Jude Children\'s Research Hospital, 2025')

        options.add_argument("assign", nargs="?")
        options.add_argument('--in','-i',action='store',required=True,nargs="*",
                            metavar='query',dest='query',
                            help='Input (multi-)fasta files to type (each contig is typed separately)')
        options.add_argument('--db','-d',action='store',required=True,nargs="*",
                            metavar='refdb',dest='refdb',
                            help='Specify viral database for the lineage classification (defaul=lineages.tsv)')
        options.add_argument('--out','-o',action='store',required=False,nargs=1,
                            metavar='outfile',dest='outfile',default="lineages.tsv",
                            help='Output file containing a summary of the assigned lineages')
        options.add_argument('--seq','-s',action='store_true',default=False,
                            dest='showseq',help='Show nucleotide sequence in the output')
        options.add_argument('--force','-f',action='store_true',default=False,
                            dest='force',help='Force overwrite output file')
        options.add_argument('--threads','-t',action='store',required=False,nargs=1,
                            metavar='threads',dest='threads',default=10,
                            help='Number of threads (default=1)')
        options.add_argument('--quiet','-q',action='store_false',default=True,
                            dest='verbose',help='Show viral lineage assignment progress')

        options=options.parse_args(args=None if sys.argv[2:] else ['--help'])

        if len(options.refdb[0:])>1:
            print("\nMultiple databases specified; only one is required!\n")
            sys.exit()
        else:
            pass

        cmdValues = {'inputSeqFiles': options.query[0:][0:],
                     'refDatabase': options.refdb[0:][0],
                     'outputLineageFile': options.outfile[0:][0] if isinstance(options.outfile,list) else options.outfile,
                     'includeSequence': options.showseq,
                     'threads': int(options.threads[0:][0]) if isinstance(options.threads,list) else options.threads,
                     'forceOverwrite': options.force, 
                     'verboseOutput': options.verbose}

        global refdbName
        refdbName = options.refdb[0:][0]

        checkDependencies()
        checkDatabase()

        return(cmdValues)

    elif sys.argv[1]=="database":
        options=argparse.ArgumentParser(sys.argv[0],
                    usage=argparse.SUPPRESS,
                    description='virallc: A tool for rapid assignment of virus lineages for given nomenclature',
                    prefix_chars='-',
                    add_help=True,
                    epilog='Written by Chrispin Chaguza, St Jude Children\'s Research Hospital, 2025')

        options.add_argument("database", nargs="?")
        options.add_argument('--showdb','-p',action='store_true',default=False,
                            dest='showdb',help='Print list of implemented databases')
        options.add_argument('--setupdb','-s',action='store_true',default=False,
                            dest='setupdb',help='Setup implemented databases')
        options.add_argument('--updatedb','-u',action='store_true',default=False,
                            dest='updatedb',help='Update implemented databases')
        options.add_argument('--version','-v',action='store_true',default=False,
                            dest='version',help='Show database version')

        if len(sys.argv)<=2:
            options=options.parse_args(args=['--help'])
        else:
            options=options.parse_args()

        if options.showdb:
            showDatabases()
            sys.exit()
        elif options.setupdb or options.updatedb:
            updateDatabases()
            sys.exit()
        elif options.version:
            showDatabaseVersions()
            sys.exit()
        else:
            sys.exit()

    elif sys.argv[1]=="version":
        global version
        print(f"\nVersion: {version}\n")
        sys.exit()

    elif sys.argv[1]=="citation":
        print("Citation information is available at https://github.com/ChrispinChaguza/virallc")
        sys.exit()

    else:
        sys.exit()

def getDatabaseName():
    global refdbName

    return(refdbName)

def getdbLocation():
    global dbpath

    return(dbpath)

def checkDependencies():
    if (not shutil.which("blastn")) or (not shutil.which("blastn")):
        print("–Install NCBI BLAST before running the program")
        sys.exit()
    else:
        pass

    if (not shutil.which("nextclade")):
        print("–Install nextclade sequence aligner before running the program")
        sys.exit()
    else:
        pass

    if (not shutil.which("mafft")):
        print("–Install mafft sequence aligner before running the program")
        sys.exit()
    else:
        pass

def updateDatabases():
    dbpath = getdbLocation()

    if os.path.exists(dbpath):
        shutil.rmtree(dbpath)
        os.makedirs(dbpath,exist_ok=True)
    else:
        os.makedirs(dbpath,exist_ok=True)

    try:
        result = os.system(f"rm -rf {dbpath} >/dev/null 2>&1")
        result = os.system(f"gitdir https://github.com/ChrispinChaguza/virallc/tree/main/viraldb >/dev/null 2>&1")
        result = os.system(f"mv viraldb {dbpath} >/dev/null 2>&1")
        result = os.system(f"rm -rf {os.path.join(dbpath,"viraldb")} >/dev/null 2>&1")

        #gitfile = fsspec.filesystem("github",org="chrispinchaguza",repo="virallc")
        #gitfile.get(gitfile.ls("viraldb"), dbpath)
    except:
        print("Failed to update databases (check the internet?)")
        sys.exit()

    sys.exit()

def showDatabases():
    dbpath = getdbLocation()
    
    if (not os.path.exists(dbpath)):
        print("Run viraL database --setupdb to setup the databases")
    else:
        for dbnum,db in enumerate(os.listdir(dbpath)):
            if os.path.isdir(os.path.join(dbpath,db)):
                if (os.path.exists(os.path.join(dbpath,db,"version.txt"))) and (os.path.exists(os.path.join(dbpath,db,"strain.db.info.txt"))):
                    dbDescription = [i.strip() for i in open(os.path.join(dbpath,db,"strain.db.info.txt"),"r")][0]
                    dbVersion = [i.strip() for i in open(os.path.join(dbpath,db,"version.txt"),"r")][0]
                    print(f"Database: {dbnum+1}; dbname: {db}; description: {dbDescription}; version: {dbVersion}; location: {os.path.join(dbpath,db)}")
                else:
                    print(f"Database: {dbnum+1}; dbname: {db}; description: None (update databases!); version: None (update databases!); location: {os.path.join(dbpath,db)}")
            else:
                pass

    sys.exit()

def showDatabaseVersions():
    dbpath = getdbLocation()
    
    if (not os.path.exists(dbpath)):
        print("Run viraL database --setupdb to setup the databases")
    else:
        for dbnum,db in enumerate(os.listdir(dbpath)):
            if os.path.isdir(os.path.join(dbpath,db)):
                if (os.path.exists(os.path.join(dbpath,db,"version.txt"))):
                    dbVersion = [i.strip() for i in open(os.path.join(dbpath,db,"version.txt"),"r")][0]
                    print(f"Database: {dbnum+1}; dbname: {db}; version: {dbVersion}")
                else:
                    print(f"Database: {dbnum+1}; dbname: {db}; version: None (update databases!)")
            else:
                pass

    sys.exit()

def getDatabaseDescription():
    dbpath = getdbLocation()
    refdbName = getDatabaseName()

    if (not os.path.exists(dbpath)):
        print("Run viraL database --setupdb to setup the databases")
    else:
        if os.path.isdir(dbpath):
            if (os.path.exists(os.path.join(dbpath,refdbName,"strain.db.info.txt"))):
                dbName = [i.strip() for i in open(os.path.join(dbpath,refdbName,"strain.db.info.txt"),"r")][0]
            else:
                print(f"File {os.path.join(dbpath,refdbName,'strain.db.info.txt')} does not exist (update databases!)")
        else:
            pass

    return(dbName)

def checkDatabase():
    refdbName = getDatabaseName()
    dbpath = getdbLocation()
   
    if (not os.path.exists(dbpath)):
        print("Run 'viraL database --setupdb' to setup the databases for lineage classification")
    else:
        if (os.path.exists(os.path.join(dbpath,refdbName))):
            if (not os.path.exists(os.path.join(dbpath,refdbName,"strain.db.info.txt"))):
                print(f"Database file {os.path.join(dbpath,refdbName,'strain.db.info.txt')} not found (update database!)")
                sys.exit()
            else:
                if (not os.path.exists(os.path.join(dbpath,refdbName,"references"))):
                    print(f"Database directory {os.path.join(dbpath,refdbName,'references')} not found (update database!)")
                    sys.exit()
                else:
                    if (not os.path.exists(os.path.join(dbpath,refdbName,"lineages"))):
                        print(f"Database directory {os.path.join(dbpath,refdbName,'lineages')} not found (update database!)")
                        sys.exit()
                    else:
                        if (not os.path.exists(os.path.join(dbpath,refdbName,"fasta"))):
                            print(f"Database directory {os.path.join(dbpath,refdbName,'fasta')} not found (update database!)")
                            sys.exit()
                        else:
                            pass
        else:
            print(f"Database directory {os.path.join(dbpath,refdbName)} not found (update database!)")
            sys.exit()

def getUniqueRandomString():
    global tmpdirpath

    randStr = os.path.join(tmpdirpath,
                f"tmpR.{random.randrange(1,1000000)}.{datetime.datetime.now()}".replace(':','').replace(' ','.').replace(',','.'))

    return(randStr)

def LoadSeqAlignment(refSeqName):
    refAlignObj = {}

    for seqItem in SeqIO.parse(refSeqName,"fasta"):
        refAlignObj[seqItem.id] = str(seqItem.seq).upper()

    return(refAlignObj)

def NumSNPsFromSeqToOthers(seqAlign,seqName):
    seqObj = {}

    if seqName in seqAlign.keys():
        seq1 = seqAlign[seqName].upper()
    else:
        return(seqObj)

    for seqItem in seqAlign.keys():
        if seqItem == seqName:
            continue

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

        seqObj[seqItem] = {"refSeqName": seqName, 
                             "seqCoverage": seqCoverage/len(seq1)*100,
                             "seqMatches": seqMatches/seqCoverage*100 if seqCoverage!=0 else 0}

    return(seqObj)

def NumSNPsBetweenSeqs(seqAlign,seqName1,seqName2):
    seqObj = {}
    seq1 = seqAlign[seqName1].upper()
    seq2 = seqAlign[seqName2].upper()
    seqCoverage = 0
    seqMatches = 0

    for nucSeq1,nucSeq2 in zip(seq1,seq2):
        if (nucSeq1 not in ["n","N","-","X"]) and (nucSeq2 not in ["n","N","-","X"]):
            seqCoverage = seqCoverage + 1

            if nucSeq1 == nucSeq2:
                seqMatches = seqMatches + 1
            else:
                continue
        else:
            continue

    seqObj = {"seqName1": seqName1, 
                 "seqName2": seqName2, 
                 "seqCoverage": seqCoverage/len(seq1)*100, 
                 "seqMatches": seqMatches/seqCoverage*100}

    return(seqObj)

def GetClosestSeq(seqAlignName,refSeqName):
    snpDistObj = NumSNPsFromSeqToOthers(seqAlignName,refSeqName)
    seqMinDist = {}

    for seqNum,seqName in enumerate(snpDistObj.keys()):
        if seqNum == 0:
            seqMinDist = {"closestSeqName": seqName, 
                             'refSeqName': refSeqName, 
                              'seqMatches': snpDistObj[seqName]['seqMatches'], 
                              'seqCoverage': snpDistObj[seqName]['seqCoverage']}
        else:
            if snpDistObj[seqName]['seqMatches'] > seqMinDist['seqMatches'] and seqName != refSeqName:
                seqMinDist = {"closestSeqName": seqName, 
                                 'refSeqName': refSeqName, 
                                 'seqMatches': snpDistObj[seqName]['seqMatches'], 
                                 'seqCoverage': snpDistObj[seqName]['seqCoverage']}
            else:
                continue

    return(seqMinDist)

def AddSeqToAlignment(seqfastafile,seqFile,refSeqFile):
    threads = 5
    outfastafile = str(getUniqueRandomString()+".aln")
    outfastafile1 = str(getUniqueRandomString()+".aln")

    nextcladecmd = f"nextclade run --silent --kmer-distance 15 --kmer-length 31 --allowed-mismatches 50 \
                       --window-size 100 --min-match-length 50 --max-alignment-attempts 5 --min-length 50 \
                       --output-fasta {outfastafile1} --input-ref {refSeqFile} {seqFile} {seqfastafile}"
    subprocess.call(nextcladecmd,shell=True,stdout=open(outfastafile1,'w'),stderr=subprocess.STDOUT)

    tmpAlignSeqs = {}

    seqName = [i.id for i in SeqIO.parse(seqFile,"fasta")][0]

    if os.path.exists(outfastafile1):
        if os.path.getsize(outfastafile1)==0:
            mafftcmd = f"mafft --quiet --auto --thread {threads} --addfull {seqFile} --keeplength {refSeqFile}"
            subprocess.call(mafftcmd,shell=True,stdout=open(outfastafile,'w'),stderr=subprocess.STDOUT)

            mafftcmd1 = f"mafft --quiet --auto --thread {threads} --addfull {seqfastafile} --keeplength {outfastafile}"
            subprocess.call(mafftcmd1,shell=True,stdout=open(outfastafile1,'w'),stderr=subprocess.STDOUT)

            os.remove(outfastafile)
        else:
            pass
    else:
        mafftcmd = f"mafft --quiet --auto --thread {threads} --addfull {seqFile} --keeplength {refSeqFile}"
        subprocess.call(mafftcmd,shell=True,stdout=open(outfastafile,'w'),stderr=subprocess.STDOUT)

        mafftcmd1 = f"mafft --quiet --auto --thread {threads} --addfull {seqfastafile} --keeplength {outfastafile}"
        subprocess.call(mafftcmd1,shell=True,stdout=open(outfastafile1,'w'),stderr=subprocess.STDOUT)

        os.remove(outfastafile)

    if os.path.exists(outfastafile1):
        tmpAlignSeqs = LoadSeqAlignment(outfastafile1)
        
        if (len(tmpAlignSeqs.keys())>0) and (seqName in tmpAlignSeqs.keys()):
            pass
        else:
            mafftcmd = f"mafft --quiet --auto --thread {threads} --addfull {seqFile} --keeplength {refSeqFile}"
            subprocess.call(mafftcmd,shell=True,stdout=open(outfastafile,'w'),stderr=subprocess.STDOUT)

            mafftcmd1 = f"mafft --quiet --auto --thread {threads} --addfull {seqfastafile} --keeplength {outfastafile}"
            subprocess.call(mafftcmd1,shell=True,stdout=open(outfastafile1,'w'),stderr=subprocess.STDOUT)

            os.remove(outfastafile)

            if os.path.exists(outfastafile1):
                if os.path.getsize(outfastafile1)==0:
                    tmpAlignSeqs = {}
                else:
                    tmpAlignSeqs1 = LoadSeqAlignment(outfastafile1)
                    if (len(tmpAlignSeqs1.keys())>0) and (seqName in tmpAlignSeqs1.keys()):
                        tmpAlignSeqs = tmpAlignSeqs1
                    else:
                        tmpAlignSeqs = {}
            else:
                tmpAlignSeqs = {}
    else:
         tmpAlignSeqs = {}

    os.remove(outfastafile1)

    return(tmpAlignSeqs)

def WriteAlignment(seqObjDict,outputfastafile):
    with open(outputfastafile,"w") as fastafile:
        for seqName in seqObjDict.keys():
            fastafile.write(">{seqName}\n{seq}\n".format(seqName=seqName,seq=seqObjDict[seqName]))

def CompareSeqsBLAST(refSeqFile,seqFile):
    tmpBlastOutFile = str(getUniqueRandomString()+".m8")
    
    blastcmd = f"blastn -query {refSeqFile} -subject {seqFile} -outfmt \
            '6 qseqid sseqid qlen slen qstart qend length pident qcovs evalue qseq sseq'"

    subprocess.call(blastcmd,shell=True,stdout=open(tmpBlastOutFile,'w'),stderr=subprocess.STDOUT)

    blastMatch = {}

    if os.path.exists(tmpBlastOutFile):
        tmpFile = [str(record).strip() for record in open(tmpBlastOutFile,"r")]
    else:
        tmpFile = []

    if len(tmpFile) == 0:
        blastcmd = f"tblastx -query {refSeqFile} -subject {seqFile} -outfmt \
                '6 qseqid sseqid qlen slen qstart qend length pident qcovs evalue qseq sseq'"
        subprocess.call(blastcmd,shell=True,stdout=open(tmpBlastOutFile,'w'),stderr=subprocess.STDOUT)

        tmpFileRec = [str(record).strip() for record in open(tmpBlastOutFile,"r")]

        if len(tmpFileRec) == 0:
            return(blastMatch)
        else:
            pass
    else:
        pass

    tmpBlastMatchtmp = [str(record).strip().split("\t") for record in open(tmpBlastOutFile,"r")]

    for i,j in enumerate(tmpBlastMatchtmp):
        tmpBlastMatchtmp[i][6] = float(tmpBlastMatchtmp[i][6])
        tmpBlastMatchtmp[i][7] = float(tmpBlastMatchtmp[i][7])
        tmpBlastMatchtmp[i][8] = float(tmpBlastMatchtmp[i][8])
        tmpBlastMatchtmp[i][9] = float(tmpBlastMatchtmp[i][9])

    tmpBlastMatch = sorted(tmpBlastMatchtmp,key=itemgetter(9), reverse=False)[0]

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

    return(blastMatch)

def GetReferenceSequence(refName):
    refdbName = getDatabaseName()
    dbpath = getdbLocation()
    
    fastafiles = {}

    for fastafile in os.listdir(os.path.join(dbpath,refdbName,"fasta")):
        if os.path.splitext(os.path.basename(fastafile))[1]==".fasta":
            fastafiles[os.path.splitext(os.path.basename(fastafile))[0]]=os.path.join(dbpath,refdbName,"fasta",fastafile)
        else:
            continue

    if refName in [i for i in fastafiles.keys()]:
        return(fastafiles[refName])
    else:
        return(False)

def GetSeqSeqForSegments(refName):
    refdbName = getDatabaseName()
    dbpath = getdbLocation()
 
    refseqFiles = {}

    for seqfile in os.listdir(os.path.join(dbpath,refdbName,"references")):
        if os.path.splitext(os.path.basename(seqfile))[1]==".fasta":
            refseqFiles[os.path.splitext(os.path.basename(seqfile))[0]]=os.path.join(dbpath,refdbName,"references",seqfile)
        else:
            continue

    if refName in [i for i in refseqFiles.keys()]:
        return(refseqFiles[refName])
    else:
        print(f"Reference sequence for {refName} not found (update database!)")
        return(False)

def GetLineageThresholds(segmentName):
    dbpath = getdbLocation()
    refdbName = getDatabaseName()

    lineageThresholds = {}
    for seqItem in [i.strip().split("\t") for i in open(os.path.join(dbpath,refdbName,"lineages","lineage.thresholds.txt"),"r")]:
        lineageThresholds[seqItem[0]] = [seqItem[1],seqItem[2]]

    return(lineageThresholds[segmentName])

def GetRefLineageNames(refName):
    dbpath = getdbLocation()
    refdbName = getDatabaseName()

    lineageFiles = {}

    for lineagefile in os.listdir(os.path.join(dbpath,refdbName,"lineages")):
        if os.path.splitext(os.path.basename(lineagefile))[1]==".tsv":
            lineageFiles[os.path.splitext(os.path.basename(lineagefile))[0]]=os.path.join(dbpath,refdbName,"lineages",lineagefile)
        else:
            continue

    if refName in [i for i in lineageFiles.keys()]:
        return(lineageFiles[refName])
    else:
        return(False)

def getReferenceMFA():
    refseqMFA = os.path.join(dbpath,refdbName,"references","references.mfa")

    return(refseqMFA)

def GetRefSeqsSegmentsForBLAST():
    refdbName = getDatabaseName()
    dbpath = getdbLocation()
    refseqs = getReferenceMFA()

    return(refseqs)

def assignLineageToSequence(segmentName,seqName):
    refLineageFile = GetRefLineageNames(segmentName)
    refLineageDB = {}

    if refLineageFile != False:
        for i in [str(i).strip().split("\t") for i in open(refLineageFile,"r")]:
            refLineageDB[i[0]] = i[1]
        if seqName in list(refLineageDB.keys()):
            return(refLineageDB[seqName])
        else:
            return(False)
    else:
        return(False)

def GetLineage(seqFastaFile):
    assignedLineageReport = {}
    blastCoverageThresholdSegment = 20

    for seqFastaItem in SeqIO.parse(seqFastaFile,"fasta"):
        tmpSeqName = seqFastaItem.id

        seqFile = str(getUniqueRandomString()+".fna")

        checkVar = False
        with open(seqFile,"w") as tmpSeqHandle:
            tmpSeqHandle.write(f">{seqFastaItem.id}____XX___\n{seqFastaItem.seq}\n")

        if len(seqFastaItem.seq)==0:
            blastResults = {}
        else:
            blastResults = CompareSeqsBLAST(GetRefSeqsSegmentsForBLAST(),seqFile)

        closestMatchedRefSeq = {}

        if len(blastResults.keys())==0 or (float(blastResults['qcovs']) < blastCoverageThresholdSegment):
            assignedLineageReport[tmpSeqName] = {"seqFastaFile": seqFastaFile, 
                                                    "seqName": tmpSeqName, 
                                                    "segment": 'NA', 
                                                    "assignedLineageName": 'NA',
                                                    "assignedLineageType": 'NA', 
                                                    "closestMatchedRefSeq": 'NA',
                                                    "closestMatchIdentity": 'NA', 
                                                    "closestMatchCoverage": 'NA',
                                                    "coverageFlag": 'No match found', 
                                                    "assignedLineageInstr": 'NA', 
                                                    "contactDetails": 'NA', 
                                                    "assignedLineageInfo": 'Insufficient query sequence coverage', 
                                                    "seq": str(seqFastaItem.seq)}
        else:
            refSegmentSeq = GetSeqSeqForSegments(blastResults['qseqid'])
            refSegmentALN = GetReferenceSequence(blastResults['qseqid'])

            outSeqAlign = AddSeqToAlignment(refSegmentALN,seqFile,refSegmentSeq)

            if len(outSeqAlign.keys())>0:
                closestMatchedRefSeq = GetClosestSeq(outSeqAlign,blastResults['sseqid'])
            else:
                pass

            assignedLineageName = ""
            assignedLineageType = ""
            assignedLineageInfo = ""

            coverageFlag = ""
            contactDetails = "chrispin.chaguza@gmail.com"

            if len(closestMatchedRefSeq.keys())==0:
                assignedLineageReport[tmpSeqName] = {"seqFastaFile": seqFastaFile,
                    "seqName": tmpSeqName, 
                    "segment": 'NA', 
                    "assignedLineageName": 'NA',
                    "assignedLineageType": 'NA', 
                    "closestMatchedRefSeq": 'NA',
                    "closestMatchIdentity": 'NA', 
                    "closestMatchCoverage": 'NA',
                    "coverageFlag": 'NA', 
                    "assignedLineageInstr": 'NA',
                    "assignedLineageInfo": "Lineage assignment failed!",
                    "contactDetails": 'NA', 
                    "seq": str(seqFastaItem.seq)}
            else:
                LineageName = assignLineageToSequence(blastResults['qseqid'],closestMatchedRefSeq['closestSeqName'])

                if float(closestMatchedRefSeq['seqMatches']) >= float(GetLineageThresholds(blastResults['qseqid'])[0]):
                    if float(closestMatchedRefSeq['seqMatches']) >= float(GetLineageThresholds(blastResults['qseqid'])[1]):
                        assignedLineageName = LineageName 
                        assignedLineageType = "Known lineage"
                        assignedLineageInfo = ""
                        assignedClosestLineage = LineageName
                    else:
                        assignedLineageName = LineageName
                        assignedLineageType = "Novel lineage?"
                        assignedClosestLineage = LineageName
                        assignedLineageInfo = "Verify with the curators"
                else:
                    assignedLineageName = LineageName
                    assignedLineageType = "Novel lineage?"
                    assignedClosestLineage = f"{LineageName} (closest)"
                    assignedLineageInfo = "Verify with the curators"

                if closestMatchedRefSeq['seqCoverage'] >= 85:
                    coverageFlag = "Excellent"
                elif closestMatchedRefSeq['seqCoverage'] >= 75:
                    coverageFlag = "Good"
                elif closestMatchedRefSeq['seqCoverage'] >= 60:
                    coverageFlag = "Fair"
                else:
                    coverageFlag = "Poor!"

                assignedLineageReport[tmpSeqName] = {"seqFastaFile": seqFastaFile, 
                    "seqName": tmpSeqName, "segment": blastResults['qseqid'], 
                    "assignedLineageName": assignedLineageName if assignedLineageName else "NA",
                    "assignedLineageType": assignedLineageType, 
                    "closestMatchedRefSeq": closestMatchedRefSeq['closestSeqName'],
                    "closestMatchIdentity": closestMatchedRefSeq['seqMatches'], 
                    "closestMatchCoverage": closestMatchedRefSeq['seqCoverage'],
                    "coverageFlag": coverageFlag, 
                    "assignedLineageInfo": assignedLineageInfo, 
                    "contactDetails": contactDetails, 
                    "seq": str(seqFastaItem.seq).upper()}

        if os.path.exists(seqFile):
            os.remove(seqFile)
        else:
            pass

    return(assignedLineageReport)

def lineageReport(seqFileName,tmpSeqFile,eachSeq,includeSequence,verbose):
    if os.path.exists(tmpSeqFile):
        assignedLineages = GetLineage(tmpSeqFile) 

        for eachSeqLineage in assignedLineages.keys():
            if includeSequence:
                lineageReport = assignedLineages[eachSeqLineage]
                lineageReportSummary = str(str("\t").join(["{seqFastaFile}",\
                                           "{seqName}",\
                                           "{segment}",\
                                           "{VirusName}",\
                                           "{assignedLineageName}",\
                                           "{assignedLineageType}",\
                                           "{closestMatchedRefSeq}",\
                                           "{closestMatchIdentity}",\
                                           "{closestMatchCoverage}",\
                                           "{coverageFlag}",\
                                           "{assignedLineageInfo}",\
                                           "{nucSequence}","\n"])).format(seqFastaFile=seqFileName,
                      seqName=lineageReport['seqName'],
                      VirusName=getDatabaseDescription(),
                      segment=lineageReport['segment'],
                      assignedLineageName=lineageReport['assignedLineageName'],
                      assignedLineageType=lineageReport['assignedLineageType'],
                      closestMatchedRefSeq=lineageReport['closestMatchedRefSeq'],
                      closestMatchIdentity=round(lineageReport['closestMatchIdentity'],2) if not \
                                              isinstance(lineageReport['closestMatchIdentity'], str) else lineageReport['closestMatchIdentity'],
                      closestMatchCoverage=round(lineageReport['closestMatchCoverage'],2) if not \
                                              isinstance(lineageReport['closestMatchCoverage'], str) else lineageReport['closestMatchCoverage'],
                      coverageFlag=lineageReport['coverageFlag'],
                      assignedLineageInfo=lineageReport['assignedLineageInfo'],
                      contactDetails=lineageReport['contactDetails'],
                      nucSequence=lineageReport['seq'])
            else:
                lineageReport = assignedLineages[eachSeqLineage]
                lineageReportSummary = str(str("\t").join(["{seqFastaFile}",\
                                           "{seqName}",\
                                           "{VirusName}",\
                                           "{segment}",\
                                           "{assignedLineageName}",\
                                           "{assignedLineageType}",\
                                           "{closestMatchedRefSeq}",\
                                           "{closestMatchIdentity}",\
                                           "{closestMatchCoverage}",\
                                           "{coverageFlag}",\
                                           "{assignedLineageInfo}","\n"])).format(seqFastaFile=seqFileName,
                      seqName=lineageReport['seqName'],
                      VirusName=getDatabaseDescription(),
                      segment=lineageReport['segment'],
                      assignedLineageName=lineageReport['assignedLineageName'],
                      assignedLineageType=lineageReport['assignedLineageType'],
                      closestMatchedRefSeq=lineageReport['closestMatchedRefSeq'],
                      closestMatchIdentity=round(lineageReport['closestMatchIdentity'],2) if not \
                                               isinstance(lineageReport['closestMatchIdentity'], str) else lineageReport['closestMatchIdentity'],
                      closestMatchCoverage=round(lineageReport['closestMatchCoverage'],2) if not \
                                               isinstance(lineageReport['closestMatchCoverage'], str) else lineageReport['closestMatchCoverage'],
                      coverageFlag=lineageReport['coverageFlag'],
                      assignedLineageInfo=lineageReport['assignedLineageInfo'],
                      contactDetails=lineageReport['contactDetails'],
                      nucSequence=lineageReport['seq'])

        if verbose:
            print(lineageReportSummary.strip())
        else:
            pass

        if os.path.exists(tmpSeqFile):
            os.remove(tmpSeqFile)
        else:
            pass

        return(lineageReportSummary)
    else:
        if os.path.exists(tmpSeqFile):
            os.remove(tmpSeqFile)
        else:
            pass

        return(f"Sequence {eachSeq.id} not processed")

def lineages():
    cmdOptionsVals = ReadCmdOptions()

    if os.path.exists(cmdOptionsVals['outputLineageFile']):
        if cmdOptionsVals['forceOverwrite']:
            os.remove(cmdOptionsVals['outputLineageFile'])
        else:
            print(f"\nOutfile '{cmdOptionsVals['outputLineageFile']}' already exists (use --force to overwrite)\n")
            sys.exit()
    else:
        pass

    global tmpdirpath

    if os.path.exists(tmpdirpath):
        shutil.rmtree(tmpdirpath)
        os.makedirs(tmpdirpath,exist_ok=True)
    else:
        os.makedirs(tmpdirpath,exist_ok=True)

    with open(cmdOptionsVals['outputLineageFile'],"w") as outputLineageHandle:
        if cmdOptionsVals['includeSequence']:
            lineageReportHeader = str("\t").join(["sequence",\
                                        "name",\
                                        "database",\
                                        "reference",\
                                        "lineage",\
                                        "category",\
                                        "closestmatch",\
                                        "identity",\
                                        "coverage",\
                                        "quality",\
                                        "description",\
                                        "query","\n"])
        else:
            lineageReportHeader = str("\t").join(["sequence",\
                                        "name",\
                                        "database",\
                                        "reference",\
                                        "lineage",\
                                        "category",\
                                        "closestmatch",\
                                        "identity",\
                                        "coverage",\
                                        "quality",\
                                        "description","\n"])

        outputLineageHandle.write(lineageReportHeader)

        if cmdOptionsVals['verboseOutput']:
            print(lineageReportHeader.strip())
        else:
            pass

        tmpList = [] 

        for seqFileName in cmdOptionsVals['inputSeqFiles']:
            if os.path.exists(seqFileName):
                for eachSeq in SeqIO.parse(seqFileName,"fasta"):
                    tmpSeqFile = str(getUniqueRandomString()+".fas")
                    
                    with open(tmpSeqFile,"w") as tmpSeqFileHandle:
                        seqName=eachSeq.id
                        seq=str(eachSeq.seq).replace("-","")
                        tmpSeqFileHandle.write(f">{seqName}\n{seq}\n")
                    
                    tmpList.append([seqFileName,tmpSeqFile,eachSeq])
            else:
                print(f"Input sequence file '{seqFileName}' does not exist")
                sys.exit()

        with concurrent.futures.ThreadPoolExecutor(cmdOptionsVals['threads']) as executor:
            jobs = [executor.submit(lineageReport,p,q,r,cmdOptionsVals['includeSequence'],cmdOptionsVals['verboseOutput']) for p,q,r in tmpList]  
            for output in concurrent.futures.as_completed(jobs): 
                outputLineageHandle.write(output.result())
   
    if os.path.exists(tmpdirpath):
        shutil.rmtree(tmpdirpath)
    else:
        pass
