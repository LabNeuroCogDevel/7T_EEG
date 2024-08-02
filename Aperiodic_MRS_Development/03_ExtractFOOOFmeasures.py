#!/usr/bin/env python3

import scipy

# General imports
import matplotlib.pyplot as plt
from matplotlib import cm, colors, colorbar

# Import MNE, as well as the MNE sample dataset
import mne
from mne import io
from mne.datasets import sample
from mne.viz import plot_topomap

# FOOOF imports
from fooof import FOOOF
from fooof import FOOOFGroup
from fooof.sim.gen import gen_power_spectrum
from fooof.sim.utils import set_random_seed
from fooof.plts.spectra import plot_spectrum
from fooof.plts.annotate import plot_annotated_model
from fooof.bands import Bands
from fooof.analysis import get_band_peak_fg
from fooof.plts.spectra import plot_spectrum
import os.path as op

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import (make_axes_locatable, ImageGrid,
                                     inset_locator)
import os.path
from os import path
import pandas as pd


def check_nans(data, nan_policy='zero'):
    """Check an array for nan values, and replace, based on policy."""

    # Find where there are nan values in the data
    nan_inds = np.where(np.isnan(data))

    # Apply desired nan policy to data
    if nan_policy == 'zero':
        data[nan_inds] = 0
    elif nan_policy == 'mean':
        data[nan_inds] = np.nanmean(data)
    else:
        raise ValueError('Nan policy not understood.')
    return data


def reject_bad_segs(raw):
    """ This function rejects all time spans annotated as 'bad_break' and concatenates the rest"""
    raw_segs = []
    for jsegment in range(len(raw.annotations) - 1):
        if raw.annotations.description[jsegment] != 'BAD_break':  # Append all other than 'bad_break'
            raw_segs.append(
                raw.copy().crop(tmin=raw.annotations.onset[jsegment], tmax=raw.annotations.onset[jsegment + 1],
                                include_tmax=False))
    return mne.concatenate_raws(raw_segs)


import glob

# load in participant files and extract the data we want
files = glob.glob(
    '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/Individual_Files/1*_2*_eyesClosed.npz') #load in eyes closed
outpath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/Individual_Files/'

overallDF = pd.DataFrame()

for fname in files:
    loaded = np.load(fname)
    sub = fname[75:89]
    cond = fname[90:100]
    subAps = loaded['aperiodic']
    offset = subAps[:, 0]
    exponent = subAps[:, 1]
    subDF = pd.DataFrame(
        {'Subject': sub, 'Channel': raw.ch_names, 'Offset': offset, 'Exponent': exponent, 'Condition': cond})

    overallDF = pd.concat([overallDF, subDF])


files = glob.glob(
    '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/Individual_Files/1*_2*_eyesOpen.npz') #load in eyes open

for fname in files:
    loaded = np.load(fname)
    sub = fname[75:89]
    cond = fname[90:98]
    subAps = loaded['aperiodic']
    offset = subAps[:, 0]
    exponent = subAps[:, 1]
    subDF = pd.DataFrame(
        {'Subject': sub, 'Channel': raw.ch_names, 'Offset': offset, 'Exponent': exponent, 'Condition': cond})

    overallDF = pd.concat([overallDF, subDF])

saveFile = outpath + 'allSubjectsFooofMeasures_20240625.csv'
overallDF.to_csv(saveFile)

# load in participant files and extract the error and rsquared values
files = glob.glob(
    '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/Individual_Files/1*_2*_eyesOpen.npz')  #load in eyes open
outpath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Results/Aperiodic_MRS_Development/Results/Individual_Files/'

## extract error and rsquared values
overallDF = pd.DataFrame()

for fname in files:
    loaded = np.load(fname)
    sub = fname[75:89]
    cond = fname[82:90]
    subError = loaded['Error']
    subRsquared = loaded['Rsquared']
    subDF = pd.DataFrame(
        {'Subject': sub, 'Channel': raw.ch_names, 'Error': subError, 'R Squared': subRsquared, 'Condition': cond})

    overallDF = pd.concat([overallDF, subDF])

files = glob.glob(
    '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/Individual_Files/1*_2*_eyesClosed.npz')  #load in eyes closed
outpath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/Individual_Files/'

## extract error and rsquared values
for fname in files:
    loaded = np.load(fname)
    sub = fname[67:81]
    cond = fname[82:92]
    subError = loaded['Error']
    subRsquared = loaded['Rsquared']
    subDF = pd.DataFrame(
        {'Subject': sub, 'Channel': raw.ch_names, 'Error': subError, 'R Squared': subRsquared, 'Condition': cond})

    overallDF = pd.concat([overallDF, subDF])

saveFile = outpath + 'allSubjectsErrorMeasures_20230516.csv'
overallDF.to_csv(saveFile)



# to extract freq specific peak info
files = glob.glob(
    '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/Individual_Files/1*_2*_eyesOpen.npz')  # change open to closed or visa versa
outpath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/Individual_Files/'
overallDF = pd.DataFrame()
for fname in files:
    loaded = np.load(fname)
    sub = fname[67:81]
    cond = fname[82:90]  # change to [82:90] for eyes closed; [71:79] for eyes open
    subPeaks = loaded[
        'peakParam']  # (centered frequency of extracted peak (CF), power of the peak over and above the aperiodic component (PW), bandwidth of extracted peak (BW))
    CF = subPeaks[:, 0]
    power = subPeaks[:, 1]
    bandwidth = subPeaks[:, 2]
    channel = subPeaks[:, 3]
    subDF = pd.DataFrame(
        {'Subject': sub, 'Channel': channel, 'Center Frequency': CF, 'Bandwidth': bandwidth, 'Power': power,
         'Condition': cond})

    overallDF = pd.concat([overallDF, subDF])

saveFile = outpath + 'allSubjectsPeakMeasures.csv'
overallDF.to_csv(saveFile)



files = glob.glob(
    '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/Individual_Files/1*_2*_eyesClosed.npz')  # change open to closed or visa versa
outpath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/Individual_Files/'
for fname in files:
    loaded = np.load(fname)
    sub = fname[67:81]
    cond = fname[82:92]  # change to [82:90] for eyes closed; [71:79] for eyes open
    subPeaks = loaded[
        'peakParam']  # (centered frequency of extracted peak (CF), power of the peak over and above the aperiodic component (PW), bandwidth of extracted peak (BW))
    CF = subPeaks[:, 0]
    power = subPeaks[:, 1]
    bandwidth = subPeaks[:, 2]
    channel = subPeaks[:, 3]
    subDF = pd.DataFrame(
        {'Subject': sub, 'Channel': channel, 'Center Frequency': CF, 'Bandwidth': bandwidth, 'Power': power,
         'Condition': cond})

    overallDF = pd.concat([overallDF, subDF])

saveFile = outpath + 'allSubjectsPeakMeasures_20230710.csv'
overallDF.to_csv(saveFile)

##run to get fooof gamma peak info
files = glob.glob(
    'H:/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/Resting_State/AfterWhole/ICAwholeClean_homogenize/1*_2*_*_Rem_rerefwhole_ICA_icapru.set')
outpath = 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/Peaks/'

# loop to find the whole spectrum
for fname in files:
    sub = fname[92:106] + '_gammaPeaks_eyesOpen.npz'
    savefile = outpath + sub
    if path.exists(savefile):
        print(f"skipping {sub}")
        continue

    raw = mne.io.read_raw_eeglab(fname)
    raw = raw.pick_types(meg=False, eeg=True, eog=False, exclude='bads')

    # raw.plot()
    raw._data = check_nans(raw._data, nan_policy='zero')
    boundary_idx = np.where(raw.annotations.description == 'boundary')[
        0]  # find where the code has triggers labeled boundary

    raw.annotations.delete(boundary_idx)  # delete boundary triggers

    raw.notch_filter(np.arange(60, 61, 1), filter_length='auto',
                     phase='zero')  # filter out the 60Hz artifact from power line noise
    # raw.plot_psd(area_mode='range', tmax=10.0, average=False);

    annotation_idx = np.where((raw.annotations.description == '16130') | (raw.annotations.description == '15362') | (
                raw.annotations.description == '1'))[
        0]
    if annotation_idx is None or len(annotation_idx) == 0:
        print(f"ERROR: no events for {fname}!")
        continue

    closed_idx = np.where((raw.annotations.description != '16130') & (raw.annotations.description != '15362') & (
                raw.annotations.description != '1'))[0]
    # find indices where the trigger DOES NOT equal 16129, 15361, 0(eyes closed) or 16130, 15262, 1 (eyes closed)
    closedData = raw.copy()  # copy the raw file
    closedData.annotations.delete(closed_idx)  # delete all other triggers

    breaks = mne.preprocessing.annotate_break(closedData, min_break_duration=5, t_start_after_previous=.1,
                                              t_stop_before_next=.1)  # find segments of data that have no annotations
    closedData.set_annotations(
        closedData.annotations + breaks)  # add in those annotations into the original file, it will be labeled as BAD_break for the segments you dont want

    appendedRawSegments = reject_bad_segs(
        closedData)  # remove segments of data you dont want and then concatenate the remaining raw files

    boundary_idx = np.where((appendedRawSegments.annotations.description == 'BAD boundary') | (
                appendedRawSegments.annotations.description == 'EDGE boundary'))[
        0]  # find indices where trigger now says bad or edge boundary. The concatenating adds this in
    appendedRawSegments_editBoundary = appendedRawSegments.copy()  # make a copy for preservation
    appendedRawSegments_editBoundary.annotations.delete(boundary_idx)  # remove triggers labeled bad or edge boundary

    appendedRawSegments_editBoundary._data = appendedRawSegments_editBoundary._data * 10e5

    spectrum = appendedRawSegments_editBoundary.compute_psd(method='welch', fmin=0, fmax=50, tmin=0, tmax=None,
                                                            picks='all', n_fft=256, n_overlap=128,
                                                            window='hamming')  # find spectrum
    psds, freqs = spectrum.get_data(return_freqs=True)  # grab frequency values corresponding to spectrum
    spectra = spectrum._data  # grab spectra values

    # Initialize a FOOOFGroup object, with desired settings
    fg = FOOOFGroup(peak_width_limits=[0.5, 12], min_peak_height=0,
                    peak_threshold=2, aperiodic_mode='fixed', max_n_peaks=4, verbose=False)

    # Define the frequency range to fit
    freq_range = [1, 50]
    # Fit the power spectrum model across all channels
    # spectraSlice = spectra[i,:,:]
    fg.fit(freqs, psds, freq_range)
    # Check the overall results of the group fits
    # fg.plot()
    # plt.show()
    # Report: fit the model, print the resulting parameters, and plot the reconstruction
    # fm = fg.get_fooof(ind=25, regenerate=True)
    # # # # Print results and plot extracted model fit
    # fm.print_results()
    # fm.plot()
    # plt.show()
    # extract specific bands
    bands = Bands({'delta': [1, 4],
                   'theta': [4, 8],
                   'alpha': [8, 13],
                   'beta': [13, 30],
                   'gamma': [30, 50]})

    # Extract gamma peaks
    gammas = get_band_peak_fg(fg,
                              bands.gamma)  # Peak data. Each row is a peak, as [CF, PW, BW]. Each row represents an individual model from the input object.
    gammas = check_nans(gammas, nan_policy='zero')
    gammas_pw = gammas[:, 1]

    sub = fname[92:106] + '_gammaPeaks_eyesOpen'
    savefile = outpath + sub
    np.savez_compressed(savefile, gammaPeaks=gammas)

# to extract freq specific peak info
files = glob.glob(
    'H:/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/Peaks/1*_2*_gammaPeaks_eyesOpen.npz')  # change open to closed or visa versa
outpath = 'H:/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/Peaks/'
overallDF = pd.DataFrame()
for fname in files:
    loaded = np.load(fname)
    sub = fname[62:76]
    cond = fname[88:96]  # change to [88:98] for eyes closed; [88:96] for eyes open
    subPeaks = loaded[
        'gammaPeaks']  # (centered frequency of extracted peak (CF), power of the peak over and above the aperiodic component (PW), bandwidth of extracted peak (BW))
    CF = subPeaks[:, 0]
    power = subPeaks[:, 1]
    bandwidth = subPeaks[:, 2]
    subDF = pd.DataFrame(
        {'Subject': sub, 'Channel': raw.ch_names, 'Center Frequency': CF, 'Bandwidth': bandwidth, 'Power': power,
         'Condition': cond})

    overallDF = pd.concat([overallDF, subDF])

saveFile = outpath + 'allSubjectsGammaPeakMeasures_20230123.csv'
overallDF.to_csv(saveFile)

