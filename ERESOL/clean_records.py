#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 11:22:20 2020

@author: ceballos
"""
from datetime import datetime
from subprocess import check_call, STDOUT
import shutil, shlex
import os



def remove_invalid_records(infile=None, ext=1, id_list=None, colname=None, outfile=None, verbose=0):
    """
    Remove invalid records from an input FITS file according to the list of Identifications
    in the provided column (it will remove the row(s) where colname==id_list)

    Parameters
    ----------
    infile : FITS
        Input fits file with all the records
    id_list : list or 1D numpy array
        Lis of identificatiions for the records to be removed
    colname : FITS column
        Column of FITS file: the id_list values refer to this column
        The row(s) where colname==id_list will be removed
    outfile : FITS file
        Output FITS file where the id_list has been removed
    verbose: verbosity level: 0 (low); 1 (chatty)

    Returns
    -------
    None. Output file is created according to "outfile" parameter

    """

    tmpFile = "pp" + str(int(datetime.timestamp(datetime.now()))) + ".fits"

    nphs = len(id_list)

    #split selection if number of removals >20 (problems in FTOOLS)
    if nphs <= 20:
        expr = "'"
        iph = 0
        while (iph < nphs-1):
            expr = expr + colname + "!= " + str(id_list[iph]) + " && "
            iph += 1
        expr = expr + colname + "!= " + str(id_list[iph]) + "'"
        comm = ("fselect infile=" + infile + "+" + str(ext) + " outfile=" + outfile +
                " clobber=yes expr=" + expr)
        try:
            print("Selecting valid records")
            if verbose:
                print(comm)
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error Selecting valid records with command:\n", comm)
            raise
    else:
        # first iteration (requires an initial "'" and input/output files are the original/final)
        expr = "'" + colname + "!= " + str(id_list[0]) + "'"
        comm = ("fselect infile=" + infile + "+" + str(ext) + " outfile=" + outfile +
                " clobber=yes expr=" + expr)
        try:
            print("Selecting valid records")
            if verbose:
                print(comm)
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error Selecting valid records with command:\n", comm)
            raise

        # other non-selections follow: create expressions with 20 selections (ftools breaks if expression is too large)
        iph = 1
        while iph < nphs:
            expr = "'"
            iiph = 1
            while (iiph < 20 and iph < nphs-1):
                expr = expr + colname + "!= " + str(id_list[iph]) + " && "
                iiph += 1
                iph  += 1

            expr = expr + colname + "!= " + str(id_list[iph]) + "'"
            iph += 1
            comm = ("fselect infile=" + outfile + "+" + str(ext) + " outfile=" + tmpFile +
                    " clobber=yes expr=" + expr)
            try:
                print("Selecting valid records")
                print("iph=" + str(iph) + "/" + str(nphs))
                if verbose:
                    print(comm)
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
            except:
                print("Error Selecting valid records with command:\n", comm)
                raise
            shutil.copy(tmpFile, outfile)
            os.remove(tmpFile)
