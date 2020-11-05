#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 12:57:00 2020

@author: ceballos

This script reads a spectral file and retrieve number of counts using pyXspec as interface

"""


import xspec
import sys
import argparse


def spec_counts(spec_file=None, emin=0., emax=0., outcts=None):
    """
    Extract spectrum counts ("net", "variance of net" or "total") using pyXspec interface

    Parameters
    ----------
    spec_file : str, optional
        Xspec spectral file. The default is None.
    emin : TYPE, optional
        Minimum energy in range (keV). The default is 0..
    emax : TYPE, optional
        Maximum energy in range (keV). The default is 0..
    outcts : TYPE, optional
        Output information about counts. The default is None.

    Returns
    -------
    Xspec counts

    """

    # read spectrum into XSPEC
    xspec.AllData.clear()
    xspec.AllModels.clear()
    s1 = xspec.Spectrum(spec_file)
    str_min = "**-" + str(float(emin))
    str_max = str(float(emax)) + "-**"
    print("Ignoring: ",str_min, str_max)
    s1.ignore(str_min)
    s1.ignore(str_max)
    (net_rate, net_rate_var, total_rate, model_rate) = s1.rate
    expos = s1.exposure
    net_cts = str(int(net_rate * expos))
    net_cts_var = str(net_rate_var * expos)
    total_cts = str(int(total_rate * expos))

    print("Net counts=", net_cts)
    print("Net counts variance=", net_cts_var)
    print("Total counts=", total_cts)

    if outcts == "netcts":
        sys.stderr.write(net_cts)
    elif outcts == "netctsvar":
        sys.stderr.write(net_cts_var)
    elif outcts == "totalcts":
        sys.stderr.write(total_cts)
    sys.exit(0)


desc = ("Retrive counts from XSPEC spectrum:" +
    "  source - bgd counts (outinfo='netcts')," +
    "  net counts variance (outinfo = 'netctsvar')" +
    "  or total counts (outinfo = 'totalcts')")
parser=argparse.ArgumentParser(description=desc)
parser.add_argument("--infile",type=str, help="name of spectral file", required=True)
parser.add_argument("--emin", type=float, help="Lower limit of energy range (keV)", required=True)
parser.add_argument("--emax", type=float, help="Upper limit of energy range (keV)", required=True)
parser.add_argument("--outinfo",type=str,choices=['netcts','netctsvar','totalcts'],
                    help="Output count information", default="totalcts")

inargs=parser.parse_args()
spec_counts(spec_file=inargs.infile,
            emin=inargs.emin,
            emax=inargs.emax,
            outcts=inargs.outinfo)
