from __future__ import print_function

def simul(nSimPulses):
    print("nSimPulses=",nSimPulses)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Create GLOBAL library from pairs of pulses', prog='simulLibsGlobal')
    parser.add_argument('--nSimPulses', type=int, help='Number of Pulses in simulated files', required=True)
    args = parser.parse_args()

    simul(nSimPulses=args.nSimPulses)

