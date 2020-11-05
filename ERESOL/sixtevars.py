#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 11:27:06 2020

@author: ceballos
"""

from os import environ
import tempfile

tmpDir = tempfile.mkdtemp()

environ["PFILES"] = tmpDir + ":" + environ["PFILES"]
environ["HEADASNOQUERY"] = ""
environ["HEADASPROMPT"] = "/dev/null/"
SIXTEinst = environ["SIXTE"] + "/share/sixte/instruments/athena-xifu/"
XIFUSIMinst = environ["SIXTE"] + "/share/xifusim/instruments/"
XMLsixte = (SIXTEinst +
            "/xifu_detector_lpa_75um_AR0.5_pixoffset_mux40_pitch275um.xml")
XMLfll = (XIFUSIMinst +
          "1pix_lpa2.5a_fll.xml")
# samplerate-dependent quantities
sampfreqs = (156250., 78125, 39062.5)  # Hz
sampids = ("", "samprate2", "samprate4")
sampStrs = ("", "_samprate2", "_samprate4")
separations = ('40000', '20000', '10000')
samplesUps = (3, 2, 2)
samplesDowns = (4, 3, 3)
nSigmss = (3.5, 4.5, 4)
preBufferPulses = 1000
scaleFactor = 0