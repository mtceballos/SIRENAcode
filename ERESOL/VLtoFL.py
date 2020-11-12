#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 13:20:58 2020

@author: ceballos
"""

from shutil import copy
import shlex
from subprocess import check_call, STDOUT

def VLtoFL(inputFile, extnum, outputFile):
    """
    Transform variable length column to fixed length using stilts
    Astropy & ftools & stilts does not work with variable-length arrays,
    so a previous conversion to fixed-format is required
    """

    print("\n Converting variable length column to fixed length column")
    print("Using:", inputFile, extnum, outputFile)
    try:
        # save inputFile into outFile
        copy(inputFile, outputFile)
        # comm = ("java -jar /home/sw/stilts.jar tpipe cmd='addcol " +
        #        "TIME2 \"(TIME*2)\"' cmd='delcols TIME2' in=" + outputFile +
        #        " out=" + outputFile)
        comm = ("stilts tpipe cmd='addcol " +
                "TIME2 \"(TIME*2)\"' cmd='delcols TIME2' ofmt=fits in=" +
                outputFile + "#" + str(extnum) + " out=" + outputFile)
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        # copy header (modified by stilts) from one FITS header to another
        # existing FITS file header
        comm = ("cphead " + inputFile + "+" + str(extnum) + " " + outputFile +
                "+" + str(extnum))
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error Coverting columns: ")
        print(comm)
        raise
