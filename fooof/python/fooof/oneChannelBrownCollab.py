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


## preprocessed data
files = glob.glob(
    '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Resting_State/AfterWhole/ICAwholeClean_homogenize/1*_2*_*_Rem_rerefwhole_ICA_icapru.set')
outpath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/fooof/'

fname = files[170] #170 for child, 14 for the adult

    raw = mne.io.read_raw_eeglab(fname, preload=True)
    raw = raw.pick_types(meg=False, eeg=True, eog=False, exclude='bads')

    # raw.plot()
    raw._data = check_nans(raw._data, nan_policy='zero')
    boundary_idx = np.where(raw.annotations.description == 'boundary')[0]  # find where the code has triggers labeled boundary

    raw.annotations.delete(boundary_idx)  # delete boundary triggers

    raw.notch_filter(np.arange(60, 61, 1), filter_length='auto', phase='zero')  # filter out the 60Hz artifact from power line noise
    # raw.plot_psd(area_mode='range', tmax=10.0, average=False);

    #make sure there are annotations
    annotation_idx = np.where((raw.annotations.description == '16130') | (raw.annotations.description == '15362') | (
                raw.annotations.description == '1'))[0]
    #if there arent any annotations then throw error message and move on to next subject
    if annotation_idx is None or len(annotation_idx) == 0:
         print(f"ERROR: no events for {fname}!")
         continue

    closed_idx = np.where((raw.annotations.description != '16130') & (raw.annotations.description != '15362') & (
                raw.annotations.description != '1'))[
        0]  #  find indices where the trigger DOES NOT equal 16129, 15361, 0(eyes closed) or 16130, 15262, 1 (eyes closed)
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
    mne.export.export_raw('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/fooof/childBrownCollab.set', appendedRawSegments_editBoundary, fmt='eeglab')




## raw data
files = glob.glob(
    '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Resting_State/remarked/1*_2*_*_Rem.set')
outpath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Jones_Collab/sampleSubs_allChannels_raw/'

fname = files[173] #173 for child, 14 for the adult

raw = mne.io.read_raw_eeglab(fname, preload=True)
raw = raw.pick_types(meg=False, eeg=True, eog=False, exclude='bads')

# raw.plot()
raw._data = check_nans(raw._data, nan_policy='zero')
boundary_idx = np.where(raw.annotations.description == 'boundary')[0]  # find where the code has triggers labeled boundary

raw.annotations.delete(boundary_idx)  # delete boundary triggers

#raw.notch_filter(np.arange(60, 61, 1), filter_length='auto', phase='zero')  # filter out the 60Hz artifact from power line noise
# raw.plot_psd(area_mode='range', tmax=10.0, average=False);


closed_idx = np.where((raw.annotations.description != '16130') & (raw.annotations.description != '15362') & (
        raw.annotations.description != '1'))[
    0]  #  find indices where the trigger DOES NOT equal 16129, 15361, 0(eyes closed) or 16130, 15262, 1 (eyes closed)
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

mne.export.export_raw('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Jones_Collab/sampleSubs_allChannels_raw/childBrownCollab_allchannels.set', appendedRawSegments_editBoundary, fmt='eeglab')
