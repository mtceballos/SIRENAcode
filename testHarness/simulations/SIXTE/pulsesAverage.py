from __future__ import print_function
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
f = fits.open('/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/tessimSPA/pp12.fits')  # mono12 but columns PIXID & PH_ID removed
tbdata = f[1].data
nrows = f[1].header['NAXIS2']
TZERO2 = f[1].header['TZERO2']
nADC = len(tbdata['ADC'][0])  # all rows are equal length
f.close()

suma = np.zeros(nADC)
listil = np.zeros(nrows)
# for each sample in record...... calculate SUMM along all the records

for il in range(0, nADC):
    # print("il=", il, "for nADC=", nADC)
    # list of ADC values in sample 'il' in all records
    listil = [x[il] for x in tbdata['ADC'][:]]
    # replace those negative (over int16 limit) by (value + 2*TZER0) (automatically done by FV)
    listil = [x+2*TZERO2 if x < 0 else x for x in listil]
    suma[il] = sum(listil)/float(nrows)
    # print("suma(", il, ")=", suma[il])


# suma = sum(tbdata['ADC'])/float(nrows)
pytstartPulse1 = 999-1  # samples starting in "1"
pytstartPulse2 = 21000-1
pulse1 = suma[pytstartPulse1:pytstartPulse1+1000]
pulse2 = suma[pytstartPulse2:pytstartPulse2+1000]
pulse12 = (pulse1 + pulse2)/2

# tbdataf = tbdata.astype([('TIME', '>f8'), ('ADC', '>f4', (2,))])

plt.plot(pulse1, 'r+')
plt.plot(pulse1, 'r--')
plt.plot(pulse2, 'bx')
plt.plot(pulse2, 'b--')
plt.plot(pulse12, 'go')
plt.plot(pulse12, 'g--')
plt.show()
