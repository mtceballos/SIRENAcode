#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 Script to convert variable-length-column fits files to fixed-length-column

Created on Wed Oct 17 12:34:48 2018

@author: ceballos
"""
import argparse
import auxpy

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Transform file with variable length cols to fixed-length',
        prog='var2fix.py')

    parser.add_argument('--inFile', help="Input file")
    parser.add_argument('--extnum', help="extension number")
    parser.add_argument('--outFile', help="Output file")

    inargs = parser.parse_args()
    print("Using:", inargs.inFile, inargs.extnum, inargs.outFile)
    auxpy.VLtoFL(inputFile=inargs.inFile,
                 extnum=inargs.extnum,
                 outputFile=inargs.outFile)
