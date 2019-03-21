"""
# PAIRS simulation
#
#      |\     |\                         |\     |\                         |\     |\
#   ___| \____| \________________________| \____| \________________________| \____| \____________
#           record                             record                             record
#   <----------------------->         <----------------------->         <----------------------->
#
# python simulPairsMultiE.py
#
#  Input parameters:
#          array (SPA|LPA1|LPA2|LPA3)
#          monoEkeV1: monochromatic energy of pulses 1(in keV)
#          monoEkeV2: monochromatic energy of pulses 2(in keV)
#          ACDC: AC or DC
#          samprate: "" or "samprate2"
#          jitter:  "" or "jitter"
#
"""

#
# --------- Read input parameters and RUN simulations -----------------
#
from __future__ import print_function
import os
import tempfile
import auxpy

cwd = os.getcwd()
tmpDir = tempfile.mkdtemp()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null/"

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Simulate pairs of pulses',
                                     prog='simulPairsMultiE.py')

    parser.add_argument('--pixName',
                        help='Extension name in pixel definition FITS file\
                        (SPA*, LPA1*, LPA2*, LPA3*)')
    parser.add_argument('--monoEnergy1', help='Monochromatic energy (keV)\
                        of input simulated first pulse')
    parser.add_argument('--monoEnergy2', help='Monochromatic energy (keV)\
                        of input simulated second pulse')
    parser.add_argument('--acbias', choices=['yes', 'no'], default="yes",
                        help='Operating Current (acbias=yes for AC\
                        or acbias=no for DC)')
    parser.add_argument('--samprate', default="", choices=['', 'samprate2'],
                        help="baseline, half_baseline")
    parser.add_argument('--jitter', default="", choices=['', 'jitter'],
                        help="no jitter, jitter")
    parser.add_argument('--noise', default="", choices=['', 'nonoise'],
                        help="noisy, no_noise")
    parser.add_argument('--stoch', default="", choices=['', 'stoch'],
                        help="no_stochastic, stochastic")
    parser.add_argument('--bbfb', default="", choices=['', 'bbfb'],
                        help="non-bbfb, bbfb")
    parser.add_argument('--pulseLength', type=int, required=True,
                        help='pulse length samples')
    parser.add_argument('--nSimPulses', type=int, required=True,
                        help='Number of simulated pulses')
    parser.add_argument('--decimation', type=int, default=1,
                        help='Decimation factor for xifusim jitter simuls')
    parser.add_argument('--separations', required=True, nargs='+',
                        help='spaced list of separations between \
                        Prim & Sec pulses')
    inargs = parser.parse_args()
    assert isinstance(inargs.separations, object)

    auxpy.simulPairs(pixName=inargs.pixName,
                     monoEkeV1=inargs.monoEnergy1,
                     monoEkeV2=inargs.monoEnergy2,
                     acbias=inargs.acbias,
                     samprate=inargs.samprate,
                     jitter=inargs.jitter,
                     noise=inargs.noise,
                     bbfb=inargs.bbfb,
                     nSimPulses=inargs.nSimPulses,
                     pulseLength=inargs.pulseLength,
                     dcmt=inargs.decimation,
                     sepsStr=inargs.separations)
