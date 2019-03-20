"""
CREATE GLOBAL LIBRARY for simulated pulses for calibration

python simulLibsGlobal.py

"""

# ----IMPORT MODULES --------------
from __future__ import print_function
import auxpy
import argparse


# ----GLOBAL VARIABLES -------------
preBufferSize = 1000

# With triggerTH=20, 0.2 keV pulses trigger 1 sample late (1001 instead
#                               of 1000), but ALL of them
#                    0.5 keV : some pulses trigger 1 sample late and some
#                               pulses trigger ok (1000) => set to 50
#                    >= 1 keV:  ALL trigger OK

#
# --------- Read input parameters and PROCESS simulations -----------------
#
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description='Create GLOBAL library from simulated pulses',
            prog='simulLibsGlobal')

    parser.add_argument('--pixName', help=('Extension name in FITS pixel \
                        definition file (SPA*, LPA1*, LPA2*, LPA3*)'),
                        required=True)
    parser.add_argument('--space', required=True,
                        choices=['ADC', 'R', 'RALL', 'RNOL', 'RFITTED'],
                        help='Input Data Space  (ADC, R, RALL, RNOL, RFITTED)')
    parser.add_argument('--samprate', default="", choices=['', 'samprate2'],
                        help="baseline, half_baseline")
    parser.add_argument('--jitter', default="", choices=['', 'jitter'],
                        help="no jitter, jitter")
    parser.add_argument('--noise', default="", choices=['', 'nonoise'],
                        help="noise, nonoise")
    parser.add_argument('--bbfb', default="", choices=['', 'bbfb'],
                        help="bbfb for stochastic tessim")
    parser.add_argument('--pulseLength', type=int, required=True,
                        help='pulse length samples')
    parser.add_argument('--nsamples', type=int, required=True,
                        help='noise samples')
    parser.add_argument('--nSimPulses', type=int, required=True,
                        help='Number of Pulses in simulated files')
    parser.add_argument('--calibEnergies', nargs='*', required=True,
                        help="list of energies (keV) from calibration files")
    parser.add_argument('--tstartPulse1All', nargs='*', type=int,
                        help="list of Tstarts for pulse 1 \
                        (if zeros, perform detection)")
    parser.add_argument('--largeFilter', type=int, default=0,
                        help='Size of extra-large filter')
    parser.add_argument('--acbias', default='yes', choices=['yes', 'no'],
                        help='Operating Current (AC )(acbias=yes) or \
                        DC (acbias=no)) [default %(default)s]')
    parser.add_argument('--createLib', type=int, required=True,
                        help='Create library (1) or only the mono files (0)')
    parser.add_argument('--noiseMat', default='no', choices=['yes', 'no'],
                        help='Should the Noise Matrices HDU be created? \
                        [default %(default)s]')
    parser.add_argument('--weightMat', default='no', choices=['yes', 'no'],
                        help='Should the Weight Matrices HDU be created? \
                        [default %(default)s]')
    parser.add_argument('--decimation', type=int, default=1,
                        help='xifusim decimation factor')

    inargs = parser.parse_args()
    len1 = 0
    if inargs.tstartPulse1All is not None:
        len1 = len(inargs.tstartPulse1All)
    lenE = len(inargs.calibEnergies)
    if len1 == 0:
        inargs.tstartPulse1All = [0 for i in range(0, lenE)]

    auxpy.simulLibsGlobal(pixName=inargs.pixName, space=inargs.space,
                          samprate=inargs.samprate, jitter=inargs.jitter,
                          noise=inargs.noise,
                          bbfb=inargs.bbfb,
                          pulseLength=inargs.pulseLength,
                          largeFilter=inargs.largeFilter,
                          libEnergies=inargs.calibEnergies,
                          tstartPulse1All=inargs.tstartPulse1All,
                          nsamples=inargs.nsamples,
                          nSimPulses=inargs.nSimPulses,
                          acbias=inargs.acbias, createLib=inargs.createLib,
                          noiseMat=inargs.noiseMat,
                          weightMat=inargs.weightMat,
                          dcmt=inargs.decimation,
                          preBufferSize=preBufferSize)
