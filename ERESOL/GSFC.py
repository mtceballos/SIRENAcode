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

# Importatnt values for GSFC channels (noise/pulse/resuduals thresholds)
#channels_ids = [ str(item) for item in list(range(1,512,2))]
channels_list = [
    {"id":"1", "noise_thS":3.5},
    {"id":"3", "noise_thS":3.5},
    {"id":"5", "noise_thS":3.5},
    {"id":"7", "noise_thS":3.5},
    {"id":"9", "noise_thS":3.65},
    {"id":"11", "noise_thS":3.5},
    {"id":"13", "noise_thS":3.5},
    {"id":"15", "noise_thS":3.5},
    {"id":"17", "noise_thS":3.5},
    {"id":"19", "noise_thS":2.6},
    {"id":"21", "noise_thS":3.5},
    {"id":"23", "noise_thS":3.7},
    {"id":"25", "noise_thS":4.0},
    {"id":"27", "noise_thS":3.7},
    {"id":"29", "noise_thS":3.8},
    {"id":"31", "noise_thS":3.5},
    {"id":"33", "noise_thS":3.5},
    {"id":"35", "noise_thS":3.5},
    {"id":"37", "noise_thS":3.5},
    {"id":"39", "noise_thS":3.5},
    {"id":"41", "noise_thS":np.NaN},
    {"id":"43", "noise_thS":3.5},
    {"id":"45", "noise_thS":3.4},
    {"id":"47", "noise_thS":3.35},
    {"id":"49", "noise_thS":3.5},
    {"id":"51", "noise_thS":3.85},
    {"id":"53", "noise_thS":3.4},
    {"id":"55", "noise_thS":4.6},
    {"id":"57", "noise_thS":3.1},
    {"id":"59", "noise_thS":3.85},
    {"id":"61", "noise_thS":3.7},
    {"id":"63", "noise_thS":3.1},
    {"id":"65", "noise_thS":3.4},
    {"id":"67", "noise_thS":3.2},
    {"id":"69", "noise_thS":3.2},
    {"id":"71", "noise_thS":3.15},
    {"id":"73", "noise_thS":3.},
    {"id":"75", "noise_thS":1.56},
    {"id":"77", "noise_thS":2.75},
    {"id":"79", "noise_thS":3.15},
    {"id":"81", "noise_thS":3.35},
    {"id":"83", "noise_thS":5.1},
    {"id":"85", "noise_thS":3.35},
    {"id":"87", "noise_thS":3.},
    {"id":"89", "noise_thS":3.3},
    {"id":"91", "noise_thS":3.15},
    {"id":"93", "noise_thS":3.3},
    {"id":"95", "noise_thS":3.15},
    {"id":"97", "noise_thS":3.5},
    {"id":"99", "noise_thS":3.5},

    {"id":"101", "noise_thS":3.1},
    {"id":"103", "noise_thS":3.4},
    {"id":"105", "noise_thS":np.NaN},
    {"id":"107", "noise_thS":2.95},
    {"id":"109", "noise_thS":3.05},
    {"id":"111", "noise_thS":3.2},
    {"id":"113", "noise_thS":5.05},
    {"id":"115", "noise_thS":3.95},
    {"id":"117", "noise_thS":3.65},
    {"id":"119", "noise_thS":3.1},
    {"id":"121", "noise_thS":3.75},
    {"id":"123", "noise_thS":2.95},
    {"id":"125", "noise_thS":3.4},
    {"id":"127", "noise_thS":3.4},
    {"id":"129", "noise_thS":3.2},
    {"id":"131", "noise_thS":3.9},
    {"id":"133", "noise_thS":3.5},
    {"id":"135", "noise_thS":3.85},
    {"id":"137", "noise_thS":3.5},
    {"id":"139", "noise_thS":4.4},
    {"id":"141", "noise_thS":2.6},
    {"id":"143", "noise_thS":4.25},
    {"id":"145", "noise_thS":3.4},
    {"id":"147", "noise_thS":3.5},
    {"id":"149", "noise_thS":4.15},
    {"id":"151", "noise_thS":3.3},
    {"id":"153", "noise_thS":2.76},
    {"id":"155", "noise_thS":4.15},
    {"id":"157", "noise_thS":4.35},
    {"id":"159", "noise_thS":3.8},
    {"id":"161", "noise_thS":3.7},
    {"id":"163", "noise_thS":3.25},
    {"id":"165", "noise_thS":4.45},
    {"id":"167", "noise_thS":3.3},
    {"id":"169", "noise_thS":np.NaN},
    {"id":"171", "noise_thS":1.3},
    {"id":"173", "noise_thS":3.1},
    {"id":"175", "noise_thS":4.6},
    {"id":"177", "noise_thS":3.2},
    {"id":"179", "noise_thS":4.2},
    {"id":"181", "noise_thS":3.5},
    {"id":"183", "noise_thS":1.64},
    {"id":"185", "noise_thS":2.46},
    {"id":"187", "noise_thS":3.8},
    {"id":"189", "noise_thS":2.0},
    {"id":"191", "noise_thS":3.75},
    {"id":"193", "noise_thS":3.5},
    {"id":"195", "noise_thS":4.5},
    {"id":"197", "noise_thS":3.25},
    {"id":"199", "noise_thS":1.5},

    {"id":"201", "noise_thS":2.2},
    {"id":"203", "noise_thS":3.6},
    {"id":"205", "noise_thS":3.2},
    {"id":"207", "noise_thS":3.5},
    {"id":"209", "noise_thS":3.5},
    {"id":"211", "noise_thS":3.55},
    {"id":"213", "noise_thS":3.2},
    {"id":"215", "noise_thS":3.35},
    {"id":"217", "noise_thS":4.6},
    {"id":"219", "noise_thS":4.15},

    {"id":"221", "noise_thS":1.8},
    {"id":"223", "noise_thS":1.6},
    {"id":"225", "noise_thS":3.85},
    {"id":"227", "noise_thS":4.4},
    {"id":"229", "noise_thS":3.4},
    {"id":"231", "noise_thS":3.3},
    {"id":"233", "noise_thS":np.NaN},
    {"id":"235", "noise_thS":1.2},
    {"id":"237", "noise_thS":3.4},
    {"id":"239", "noise_thS":4.1},
    {"id":"241", "noise_thS":3.4},
    {"id":"243", "noise_thS":3.7},
    {"id":"245", "noise_thS":3.7},
    {"id":"247", "noise_thS":3.85},
    {"id":"249", "noise_thS":5.3},
    {"id":"251", "noise_thS":3.35},
    {"id":"253", "noise_thS":1.86},
    {"id":"255", "noise_thS":4.55},
    {"id":"257", "noise_thS":3.2},
    {"id":"259", "noise_thS":4.1},
    {"id":"261", "noise_thS":3.05},
    {"id":"263", "noise_thS":3.2},
    {"id":"265", "noise_thS":3.27},
    {"id":"267", "noise_thS":3.4},
    {"id":"269", "noise_thS":3.2},
    {"id":"271", "noise_thS":3.5},
    {"id":"273", "noise_thS":2.2},
    {"id":"275", "noise_thS":2.3},
    {"id":"277", "noise_thS":3.35},
    {"id":"279", "noise_thS":3.35},
    {"id":"281", "noise_thS":1.68},
    {"id":"283", "noise_thS":3.4},
    {"id":"285", "noise_thS":4.},
    {"id":"287", "noise_thS":3.5},
    {"id":"289", "noise_thS":3.7},
    {"id":"291", "noise_thS":1.74},
    {"id":"293", "noise_thS":3.85},
    {"id":"295", "noise_thS":2.9},
    {"id":"297", "noise_thS":np.NaN},
    {"id":"299", "noise_thS":3.25},

    {"id":"301", "noise_thS":3.5},
    {"id":"303", "noise_thS":4.15},
    {"id":"305", "noise_thS":3.2},
    {"id":"307", "noise_thS":1.42},
    {"id":"309", "noise_thS":3.65},
    {"id":"311", "noise_thS":4.45},
    {"id":"313", "noise_thS":1.66},
    {"id":"315", "noise_thS":3.65},
    {"id":"317", "noise_thS":4.4},
    {"id":"319", "noise_thS":3.1},
    {"id":"321", "noise_thS":3.05},
    {"id":"323", "noise_thS":3.6},
    {"id":"325", "noise_thS":3.15},
    {"id":"327", "noise_thS":3.15},
    {"id":"329", "noise_thS":2.8},
    {"id":"331", "noise_thS":3.05},
    {"id":"333", "noise_thS":3.05},
    {"id":"335", "noise_thS":3.1},
    {"id":"337", "noise_thS":3.8},
    {"id":"339", "noise_thS":2.95},
    {"id":"341", "noise_thS":2.75},
    {"id":"343", "noise_thS":3.1},
    {"id":"345", "noise_thS":3.4},
    {"id":"347", "noise_thS":1.46},
    {"id":"349", "noise_thS":1.36},
    {"id":"351", "noise_thS":1.25},
    {"id":"353", "noise_thS":3.15},
    {"id":"355", "noise_thS":3.35},
    {"id":"357", "noise_thS":5.1},
    {"id":"359", "noise_thS":3.5},
    {"id":"361", "noise_thS":np.NaN},
    {"id":"363", "noise_thS":3.2},
    {"id":"365", "noise_thS":3.1},
    {"id":"367", "noise_thS":3.07},
    {"id":"369", "noise_thS":3.2},
    {"id":"371", "noise_thS":3.15},
    {"id":"373", "noise_thS":3.15},
    {"id":"375", "noise_thS":3.2},
    {"id":"377", "noise_thS":3.75},
    {"id":"379", "noise_thS":1.66},
    {"id":"381", "noise_thS":3.28},
    {"id":"383", "noise_thS":3.7},
    {"id":"385", "noise_thS":1.32},
    {"id":"387", "noise_thS":4.15},
    {"id":"389", "noise_thS":3.6},
    {"id":"391", "noise_thS":4.65},
    {"id":"393", "noise_thS":4.35},
    {"id":"395", "noise_thS":4.45},
    {"id":"397", "noise_thS":3.5},
    {"id":"399", "noise_thS":4.05},

    {"id":"401", "noise_thS":3.5},
    {"id":"403", "noise_thS":3.85},
    {"id":"405", "noise_thS":3.3},
    {"id":"407", "noise_thS":3.25},
    {"id":"409", "noise_thS":3.65},
    {"id":"411", "noise_thS":3.75},
    {"id":"413", "noise_thS":3.95},
    {"id":"415", "noise_thS":3.95},
    {"id":"417", "noise_thS":1.5},
    {"id":"419", "noise_thS":4.2},
    {"id":"421", "noise_thS":4.4},
    {"id":"423", "noise_thS":3.4},
    {"id":"425", "noise_thS":np.NaN},
    {"id":"427", "noise_thS":3.5},
    {"id":"429", "noise_thS":3.4},
    {"id":"431", "noise_thS":4.2},
    {"id":"433", "noise_thS":4.05},
    {"id":"435", "noise_thS":3.85},
    {"id":"437", "noise_thS":3.95},
    {"id":"439", "noise_thS":4.1},
    {"id":"441", "noise_thS":2.4},
    {"id":"443", "noise_thS":4.0},
    {"id":"445", "noise_thS":3.8},
    {"id":"447", "noise_thS":3.65},
    {"id":"449", "noise_thS":5.1},
    {"id":"451", "noise_thS":4.95},
    {"id":"453", "noise_thS":2.9},
    {"id":"455", "noise_thS":3.85},
    {"id":"457", "noise_thS":3.6},
    {"id":"459", "noise_thS":1.38},
    {"id":"461", "noise_thS":4.65},
    {"id":"463", "noise_thS":3.95},
    {"id":"465", "noise_thS":4.5},
    {"id":"467", "noise_thS":4.1},
    {"id":"469", "noise_thS":3.5},
    {"id":"471", "noise_thS":3.6},
    {"id":"473", "noise_thS":3.7},
    {"id":"475", "noise_thS":4.05},
    {"id":"477", "noise_thS":3.7},
    {"id":"479", "noise_thS":3.95},
    {"id":"481", "noise_thS":3.8},
    {"id":"483", "noise_thS":4.4},
    {"id":"485", "noise_thS":4.03},
    {"id":"487", "noise_thS":3.9},
    {"id":"489", "noise_thS":np.NaN},
    {"id":"491", "noise_thS":3.95},
    {"id":"493", "noise_thS":3.95},
    {"id":"495", "noise_thS":4.25},
    {"id":"497", "noise_thS":3.9},
    {"id":"499", "noise_thS":4.05},

    {"id":"501", "noise_thS":3.95},
    {"id":"503", "noise_thS":4.0},
    {"id":"505", "noise_thS":4.3},
    {"id":"507", "noise_thS":3.35},
    {"id":"509", "noise_thS":4.45},
    {"id":"511", "noise_thS":3.9}
]
channels_noise_thS = np.array([channel_dict["noise_thS"] for channel_dict in channels_list])
channels_ids = [channel_dict["id"] for channel_dict in channels_list]
