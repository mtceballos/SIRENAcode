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
          samprate: "" or "samprate2" or "samprate4"
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
    parser.add_argument('--nSimPulses', type=int, required=True,
                        help='Number of simulated pulses')
    parser.add_argument('--pulseLength', type=int, required=True,
                        help='pulse length samples')
    parser.add_argument('--acbias', choices=['yes', 'no'], default="yes",
                        help='Operating Current (acbias=yes for AC or '
                        'acbias=no for DC)')
    parser.add_argument('--samprate', default="",
                        choices=['', 'samprate2'],
                        help="baseline, half_baseline")
    parser.add_argument('--jitter', default="",
                        choices=['', 'jitter'],
                        help="no jitter, jitter")
    parser.add_argument('--noise', default="",
                        choices=['', 'nonoise'],
                        help="noise, no_noise")
    parser.add_argument('--bbfb', default="",
                        choices=['', 'bbfb'],
                        help="non-bbfb, bbfb")
    parser.add_argument('--decimation', type=int, default=1,
                        help='Decimation factor for xifusim jitter simuls')

    inargs = parser.parse_args()

    auxpy.simulSingles(pixName=inargs.pixName,
                       monoEkeV=inargs.monoEnergy,
                       acbias=inargs.acbias,
                       samprate=inargs.samprate,
                       jitter=inargs.jitter,
                       noise=inargs.noise,
                       bbfb=inargs.bbfb,
                       nSimPulses=inargs.nSimPulses,
                       pulseLength=inargs.pulseLength,
                       dcmt=inargs.decimation)
