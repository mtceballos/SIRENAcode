#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:30:52 2020

@author: GSFC + MTC
"""
import numpy as np
import math
import matplotlib.pyplot as plt


def averagePulse(data, thresh=None, preBuffer=100, numSamples=5000):
    """
    Averages pulse records.
    Aligment is performed before averaging

    data: numpy 2D array with list of pulses (one in every row)
    thresh: (float) initial threshold to align starting points of pulses.
            If no 'thresh' is provided, pulses are not previously aligned
    nbase: (int) number of samples to calculate the mean baseline to subtract to the average pulse
    numSamples: (int) size of average pulse
    """
    npulses = data.shape[0]
    ave_pulse = np.zeros(numSamples)
    nvalid = 0
    for ip in range(npulses):
        if thresh:
            s0 = np.where(data[ip,:] > thresh)[0][0]
            stSample = s0-preBuffer
        else:
            stSample = 0
        fnSample = stSample + numSamples
        if stSample < 0:
            print("Possible pulse in preBuffer:")
            print("      For ip=",ip,"s0,ini,fin=", s0,stSample,fnSample)
            print("      Pulse not considered in average")
            continue
        nvalid += 1
        fnSample = stSample + numSamples
        ave_pulse += data[ip,stSample:fnSample]

    factor = 1./nvalid
    ave_pulse *= factor
    nbase = int(0.8 * preBuffer)
    meanBase = np.mean(ave_pulse[:nbase])
    ave_pulse -= meanBase
    return ave_pulse


def autoDeterminePulseWindowAndThreshold(records,nbase,nsigma, numSamples=2000, plot=False, ax=None):
    """
    Determnine pulse window and pulse threshold
    Parameters
    ----------
    records : numpy 2D array [irecord, record_data]
        numpy 2D array with list of pulses (one in every row)
    nbase : int
        Number of samples of prebuffer to calculate mean baseline
    nsigma : float
        Number of sigmas above mean baseline to put the final threshold
    numSamples: int
        Number of samples to define the average pulse over which calculation of window is done
    plot : bool
        Plotting of average pulse and pulse window
    ax : axis
        plot axis for an already opened figure

    Raises
    ------
    RuntimeError
        DESCRIPTION.

    Returns
    -------
    istart : int
        Initital sample of window
    iend : int
        Final sample of window
    threshold :
        Final threshold that defines pulse window.
    avg : numpy 1D array
        Average record

    """
    # Compute average record (mean baseline subtracted)
    th0 = 0.4*np.max(records[0,:])
    #print("th0=",th0)
    avg = averagePulse(data=records, thresh=th0, preBuffer=1000, numSamples=numSamples)

    # and average baseline and baseline^2.
    nrecords = records.shape[0]
    basesum1 = 0.
    basesum2 = 0.
    #avg=np.zeros(numSamples,dtype=np.float32)
    for irec in range(nrecords):
        #avg += records[irec, :numSamples]
        basesum1 += np.sum(records[irec,:nbase])
        basesum2 += np.sum(records[irec,:nbase]**2)
        #~ print('{:12.1f}'.format(basesum1),'{:12.1f}'.format(basesum2))
    factor=1./nrecords
    basesum1 *= (factor/nbase)
    basesum2 *= (factor/nbase)
    #avg*=factor
    #avg-=basesum1
    # print(basesum1,basesum1**2,basesum2,basesum2-basesum1**2)

	# Determine pulse threshold from baseline standard deviation
	# and find where average pulse crosses it.
    threshold=nsigma*math.sqrt(basesum2-basesum1**2)
    istart=0

    while (istart < numSamples) and (avg[istart] < threshold):
        istart+=1
    iend=istart+1
    while (iend < numSamples) and (avg[iend] > threshold):
        iend+=1
    if iend >= numSamples:
        raise RuntimeError("No valid pulse window found.")

    if plot:
        #fig = plt.figure(figsize=(6,3))
        #ax = fig.add_subplot(1, 1, 1)
        ax.plot(range(numSamples), avg)
        ax.axvline(istart, ls="--", color="gray")
        ax.axvline(iend, ls="--", color="gray")
        ax.axhline(threshold, color="tab:green", ls="--")
        ax.set_xlabel("record sample")
        ax.set_ylabel("ADC amplitude (a.u.)")
        ax.set_title("Average pulse and window")


    return (istart, iend, threshold, avg)

def categorize(record0, istart, iend, thresh, joff):
    """
    Catagorize record: pulse, noise, reject
    It checks samples from istart to iend to see whether they are above or below threshold

    Parameters
    ----------
    record0 : numpy 1D array
        numpy array with values of BACKGROUND SUBTRACTED record
    istart : int
        starting sample of average pulse window
    iend : int
        final sample of average pulse window
    thresh : int
       threshold level used to define the window (and to align pulses with average pulse)
    joff : int
        jitter offset so as not to check pulse values very close to istart because of jitter

    Returns
    -------
    out : dictionary
        classification of record: pulse, noise, rejected
        out={}
		out["noise"]=0
		out["pulse"]=0
		out["rejected"]=0
        out["rejected_comm"]=0

    """

	# Prepare output
    out={}
    out["noise"]=0
    out["pulse"]=0
    out["rejected"]=0
    out["rejected_comm"]=0


    # align record and average record at 'istart' sample
    s0 = np.where(record0 > thresh)[0][0]
    stSample = s0 - istart
    if stSample < 0: # record has a preBuffer shorter than istart
        iend = s0 + (iend -istart)
        istart = s0
        trec = np.copy(record0)
    else: # record has a preBuffer larger than istart
        trec = record0[stSample:]

    #trec = np.copy(record0)
	# Record is noise if below noise threshold at all samples
    if abs(np.amin(trec)) < thresh and abs(np.amax(trec)) < thresh:
        out["noise"]=1
        return out

	# Record is a pulse if:
	#  1. below the pulse threshold before istart
	#  2. rises above the pulse threshold in [istart:iend]
	#  3. does not rise above (can fall below) pulse threshold beyond iend
	#  4. does not fall below lower noise threshold at any time

    # give a margin because of jitter
    if np.amax(trec[:(istart-joff)]) > thresh:
        out["rejected"]=1  # initial undetected pulse or high excursion
        out["rejected_comm"] = ("Above threshold befor istart")
        return out
    if np.amax(trec[(istart-joff):iend+1]) < thresh:
        out["rejected"]=1  # pulse goes below threshold unexpectedly
        out["rejected_comm"] = "Max inside pulse is below threshold"
        return out
    for i in range(iend+4,len(trec)):
        if trec[i] > thresh and trec[i] > trec[i-1]:
            out["rejected"]=1 # pulse in the tail?
            out["rejected_comm"] = ("Above thres in i=" + str(i) + " and rec[i]>rec[i-1]")
            #print("For i=",i, " record=", trec[i])
            return out
    if abs(np.amin(trec[iend+1:])) > thresh:
        out["rejected"]=1  # tail goes down below threshold of lower noise
        out["rejected_comm"] = ("Min rec past iend+1 is above threshold")
        return out
    out["pulse"]=1
    return out

def chi2(data,model, std=1.):
    chisq = np.sum(((data-model)/std)**2)
    return chisq

def deviatonFromAveragePulse(record, aveRecord, sigma=None, thresh=None):
    """
    Calculate chisq2 deviation from average pulse

    Parameters
    ----------
    record : numpy 1D array
        Numpy 1D array with record to be checked
    aveRecord : numpy 1D array
        Average record
    thresh : (float)
        Threshold value to align record and average record
        Default value is None (no alignment is required)

    Returns
    -------
    Chi square value of comparison

    """

    # alignment of record and average record
    reclen = len(record)
    avelen = len(aveRecord)
    s0_rec = np.where(record > thresh)[0]
    s0_ave = np.where(aveRecord > thresh)[0]
    maxlen = min(s0_rec,s0_ave) + min((reclen-s0_rec), (avelen-s0_ave))
    if s0_rec <= s0_ave:
        new_ave = aveRecord[(s0_ave-s0_rec):maxlen]
        new_rec = np.copy(record[:maxlen])
    else:
        new_ave = np.copy(aveRecord[:maxlen])
        new_rec = record[(s0_rec-s0_ave):maxlen]

    chisq = chi2(new_rec, new_ave, sigma)
    return chisq


def nrecs_larger_than(data, val):
    """
    Calculate nunmber of rows in which data is larger than 'val'

    Parameters:
    -----------
    data: !D numpy array of input data
    val: (1D float list/np array) threshold value(s)

    Returns
    -------

    Number of instances where data>val foreach value in val
    """

    if isinstance(val, float):
        nn = len(data[data>val])
    else:
        nvals = len(np.array(val))
        nn = np.zeros(nvals)

        for i in range(nvals):
            nn[i] = len(data[data>val[i]])
    return nn
