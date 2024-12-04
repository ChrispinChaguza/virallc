import argparse
import sys
import os

def ReadCmdOptions():
    options=argparse.ArgumentParser(sys.argv[0],
            description='RotaL: A tool for assigning rotavirus lineages',
        prefix_chars='-',
        add_help=True,
        epilog='Written by Chrispin Chaguza, Yale School of Public Health, 2023')
    options.add_argument('--sequences','-s',action='store',required=True,nargs="*",
        metavar='inputSequences',dest='inputSequences',help='Input (multi-)fasta files to type (each contig is typed separately)')
    options.add_argument('--csv','-c',action='store_true',default=False,
        dest='csvOutput',help='Write output in CSV format (default: TSV)')
    options.add_argument('--output','-o',action='store',required=True,nargs=1,
        metavar='outputLineageFile',dest='outputLineageFile',help='Output file containing a summary of the assigned lineages')
    options.add_argument('--nucseq','-n',action='store_true',default=False,
        dest='includeSequence',help='Show nucleotide sequence in the output')
    options.add_argument('--updatedb','-u',action='store_true',default=False,
        dest='updatedb',help='Setup or update the rotavirus lineage classification database')
    options.add_argument('--quiet','-q',action='store_false',default=True,
        dest='verboseOutput',help='Show the lineage assignment progress')

    options=options.parse_args()

    inputSeqFiles = options.inputSequences[0:][0:]
    outputLineageFile = options.outputLineageFile[0:][0]
    csvOutput = "\t"
    if options.csvOutput:
        csvOutput = ","
    else:
        csvOutput = "\t"

    verboseOutput = options.verboseOutput
    includeSequence = options.includeSequence
    updatedb = options.updatedb

    cmdValues = {'inputSeqFiles': inputSeqFiles,'outputLineageFile': outputLineageFile,'csvOutput': csvOutput,
            'includeSequence': includeSequence, 'updatedb': updatedb,'verboseOutput': verboseOutput}

    return(cmdValues)



