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
from .SeqLib import lineages

def main():
    lineages()            

if __name__=="__main__":
    main()
