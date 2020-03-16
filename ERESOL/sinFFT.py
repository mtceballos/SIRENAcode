#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 17:34:02 2019

@author: ceballos
"""

import matplotlib.pyplot as plt
import numpy as np
#%matplotlib inline

f = 25  # Frequency, in cycles per second, or Hertz
f_s = 100  # Sampling rate, or number of measurements per second

# define a duration time stream of 2 secs
t = np.linspace(0, 2, 2 * f_s, endpoint=False)
x = np.sin(f * 2 * np.pi * t)

fig, ax = plt.subplots()
ax.plot(t, x)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Signal amplitude')
plt.show()
















# =============================================================================
# Amp = 7
# Fs = 500  # sampling rate
# Ts = 1.0/Fs  # sampling interval
# t = np.arange(0, 2, Ts)  # time vector
# 
# ff = 25   # frequency of the signal
# y = Amp*np.sin(2*np.pi*ff*t)
# 
# n = len(y)  # length of the signal
# print("Signal length n=", n)
# # k = np.arange(n)
# # T = n/Fs
# # print("T=", T)
# # frq = k/T  # two sides frequency range = range(n)
# # print(frq)
# frq = range(n)
# frqpos = frq[1:int(n/2)]  # one side frequency range
# 
# Y = np.fft.fft(y)
# Ynorm = np.fft.fft(y)/n  # fft computing and normalization
# Ypos = 2*Y[1:int(n/2)]/n
# YnegYpos = np.zeros(n)
# YnegYpos = np.append(Ynorm[int(n/2):n], np.append(0, Ynorm[1:int(n/2)]))
# Y_csd = np.sqrt(abs(YnegYpos)**2) * np.sqrt(2*Ts/n)
# 
# fig, ax = plt.subplots(3, 1)
# ax[0].plot(t, y)
# ax[0].set_xlabel('Time')
# ax[0].set_ylabel('Amplitude')
# # plotting the (positive and total) spectrum   (ok)
# ax[1].plot(frqpos, abs(Ypos), 'r', color="blue", linestyle="--")
# ax[1].plot(range(int(-n/2), int(n/2)), abs(YnegYpos), 'r')
# ax[1].set_ylabel('|Y(freq)|')
# ax[1].axvline(x=-25, color="gray", linestyle="--")
# ax[1].axvline(x=25, color="gray", linestyle="--")
# 
# ax[2].plot(range(int(-n/2), int(n/2)), Y_csd, 'r')  # plotting the spectrum
# ax[2].set_xlabel('Freq (Hz)')
# ax[2].set_ylabel('CSD')
# plt.show()
# 
# # CSD: sqrt(E(|Y(freq)|^2))
# 
# =============================================================================
