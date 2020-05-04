#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 11:17:30 2020

@author: ceballos
"""


def tsfn(i, nproc, comm, simTimeN):
    """
    Funcion to run a single process of noise generation
    """
    if i == 0:
        tstart = 0.
    else:
        tstart = i * float(simTimeN)/10. + 0.01
    tend = float(simTimeN)/10. * (1 + i)
    outfile = "tmpbbfb" + str(i) + ".fits"
    comm += (" Streamfile=" + outfile +
             " tstart=" + str(tstart) + " tstop=" + str(tend))

    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error creating ForNoise i=", i)
        raise
    updateHISTORY(outfile, comm)
    rmLastAndFirst(outfile, 0)
