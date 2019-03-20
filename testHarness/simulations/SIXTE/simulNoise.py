"""
# NOISE spectrum simulation (with xifusim)
#
# python simulNoise.py
#
#  Input parameters:
#          pixName (SPA*|LPA1*|LPA2*|LPA3*)
#          pulseLength: length of pulses (in samples)
*          nsamples: samples for the noise
#          space: ADC (current space) or R or RALL or RNOL or RFITTED
#                                        (resistance space)
#          simTimeN: simulation time (s)
#
#
#
# 1) Simulate 10s stream with no events to calculate Baseline
#                             (pixdetillum + xifusim) --> not requiered anymore
# 2) Simulate 100s stream with no events as input for gennoisespec
#                             (tesconstpileup + xifusim)
# 3) Obtain noise spectrum with gennoisespec
#
"""

# ----IMPORT MODULES --------------
from __future__ import print_function
import os
import argparse
from auxpy import simulNoise
import tempfile

Ifit = 45.3E-6
tmpDir = tempfile.mkdtemp()
tmpFile = tempfile.TemporaryFile()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]

#
# --------- Read input parameters and PROCESS simulations -----------------
#
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description='Create NOISE spectrum', prog='simulNoise')

    parser.add_argument('--pixName', required=True,
                        help='Extension name in FITS pixel definition file:\
                        SPA*, LPA1*, LPA2*, LPA3*')
    parser.add_argument('--samprate', default="", choices=['', 'samprate2'],
                        help='baseline, half_baseline')
    parser.add_argument('--jitter',  default="", choices=['', 'jitter'],
                        help='no_jitter, jitter')
    parser.add_argument('--bbfb', default="", choices=['', 'bbfb'],
                        help="dobbfb=n, dobbfb=y")
    parser.add_argument('--pulseLength', required=True, type=int,
                        help='pulse length in samples')
    parser.add_argument('--space', required=True,
                        choices=['ADC', 'R', 'RALL', 'RNOL', 'RFITTED'],
                        help='Data space: ADC for current, R or RALL or\
                        RNOL or RFITTED for resistance')
    parser.add_argument('--acbias', default="yes",
                        help='AC (acbias=yes) or DC (acbias=no)\
                        [default %(default)s]')
    parser.add_argument('--scaleFactor', default=0.0, type=float,
                        help='Param for gennoise[default %(default)s]')
    parser.add_argument('--samplesUp', default=2, type=int,
                        help='Param samplesUp for gennoise [def %(default)s]')
    parser.add_argument('--nSgms', default=5., type=float,
                        help='Param nSgms for gennoise [default %(default)s]')
    parser.add_argument('--simTimeN', default=100,
                        help='Simulation time (s) for noise spectra\
                        calculation [default %(default)s]')
    parser.add_argument('--nintervals', default=1000, type=int,
                        help='Number of intervals in gennoisespec for spectra\
                        calculation [default %(default)s]')
    parser.add_argument('--decimation', type=int, default=1,
                        help='xifusim decimation factor')

    inargs = parser.parse_args()

    simulNoise(pixName=inargs.pixName, samprate=inargs.samprate,
               jitter=inargs.jitter, bbfb=inargs.bbfb,
               pulseLength=inargs.pulseLength,
               space=inargs.space, acbias=inargs.acbias,
               scaleFactor=inargs.scaleFactor,
               samplesUp=inargs.samplesUp,
               nSgms=inargs.nSgms, nintervals=inargs.nintervals,
               simTimeN=inargs.simTimeN,
               dcmt=inargs.decimation)
