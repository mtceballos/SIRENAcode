"""
 SINGLES simulation


      |\       |\       |\       |\       |\     |\
   ___| \______| \______| \______| \______| \____| \______
    record   record   record
   <-------><-------><------->

 python simulSingles.py

  Input parameters:
          array (SPA|LPA1|LPA2|LPA3)
          monoEkeV: monochromatic energy of pulses 1(in keV)
          ACDC: AC or DC
          samprate: "" or "samprate2"
          jitter:  "" or "jitter"
          noise:   "" or "nonoise"
"""


#
# --------- Read input parameters and RUN simulations -----------------
#
import os
import tempfile
import auxpy

cwd = os.getcwd()
tmpDir = tempfile.mkdtemp()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null/"

#nSimPulses = 20000
#nSimPulses = 1000
singleSeparation = 40000  # separation from secondary-->primary for next record

pixel = 1
preBufferSize = 1000
# with these thresholds, 0.2keV pulses are triggered at 999, 0.5keV @999
# or @1000 and larger pulses @10000


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
            description='Simulate pairs of pulses', prog='simulSingles')

    parser.add_argument('--pixName',
                        help='Extension name in pixel definition FITS file \
                        (SPA*, LPA1*, LPA2*, LPA3*)')
    parser.add_argument('--monoEnergy',
                        help='Monochromatic energy (keV) of input \
                        simulated pulse') 
    parser.add_argument('--nSimPulses',type=int, required=True,
                        help='Number of simulated pulses')
    parser.add_argument('--pulseLength', type=int, required=True,
                        help='pulse length samples')
    parser.add_argument('--acbias', choices=['yes', 'no'], default="yes",
                        help='Operating Current (acbias=yes for AC or '
                        'acbias=no for DC)')
    parser.add_argument('--samprate', default="", choices=['', 'samprate2'],
                        help="baseline, half_baseline")
    parser.add_argument('--jitter', default="", choices=['', 'jitter'],
                        help="no jitter, jitter")
    parser.add_argument('--noise', default="", choices=['', 'nonoise'],
                        help="noise, no_noise")
    parser.add_argument('--stoch', default="", choices=['', 'stoch'],
                        help="non-stochastic, stochastic")
    parser.add_argument('--bbfb', default="", choices=['', 'bbfb'],
                        help="non-bbfb, bbfb")

    inargs = parser.parse_args()
    auxpy.simulSingles(pixName=inargs.pixName, monoEkeV=inargs.monoEnergy,
                       acbias=inargs.acbias, samprate=inargs.samprate,
                       jitter=inargs.jitter, noise=inargs.noise,
                       stoch=inargs.stoch, bbfb=inargs.bbfb,
                       nSimPulses=inargs.nSimPulses,
                       singleSeparation=singleSeparation,
                       pixel=pixel, preBufferSize=preBufferSize,
                       pulseLength=inargs.pulseLength)
