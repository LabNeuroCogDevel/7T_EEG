# Import some general scientific python libraries
import numpy as np
import pandas as pd

import glob
import mne
import matplotlib.pyplot as plt

from os import path
import pandas as pd

from fooof.analysis import get_band_peak_fg

# Import the parameterization model objects
from specparam import SpectralModel, SpectralGroupModel, fit_models_3d

# Import useful parameterization related utilities and plot functions
from specparam.bands import Bands
from specparam.analysis import get_band_peak_group
from specparam.utils import trim_spectrum
from specparam.utils.data import subsample_spectra
from specparam.sim.gen import gen_aperiodic
from specparam.data import ModelSettings
from specparam.plts.templates import plot_hist
from specparam.plts.spectra import plot_spectra
from specparam.plts.periodic import plot_peak_fits, plot_peak_params
from specparam.plts.aperiodic import plot_aperiodic_params, plot_aperiodic_fits
from specparam.analysis import get_band_peak
from specparam.objs import fit_models_3d, combine_model_objs

# Import simulation & IO utilities to help with the example
from specparam.sim import sim_group_power_spectra
from specparam.sim.utils import create_freqs
from specparam.sim.params import param_sampler
from specparam.objs.utils import average_group, combine_model_objs, compare_model_objs
from specparam.bands import Bands

from specparam.utils.io import load_group_model

# Import functions to examine frequency-by-frequency error of model fits
from specparam.analysis.error import compute_pointwise_error, compute_pointwise_error_group
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


# Example Code
freq_range = [3, 40]
freq_res = 0.25

# Set up the shape of the data
n_conditions = 3
n_channels = 10
n_freqs = len(create_freqs(freq_range, freq_res))

# Define parameters for the simulated power spectra
ap_opts = param_sampler([[0, 1.0], [0, 1.5], [0, 2]])
pe_opts = param_sampler([[], [10, 0.25, 1], [10, 0.25, 1, 20, 0.15, 1]])
# Simulate power spectra, and organize into a 3D array
spectra = []
for ind in range(n_conditions):
    freqs, powers = sim_group_power_spectra(n_channels, freq_range, ap_opts,
                                            pe_opts, freq_res=freq_res)
    spectra.append(powers)

# Convert collected spectra into a numpy array
spectra = np.array(spectra)

# Initialize a SpectralGroupModel object, with desired settings
fg = SpectralGroupModel(peak_width_limits=[1, 6], min_peak_height=0.1)
# Fit the 3D array of power spectra
fgs = fit_models_3d(fg, freqs, spectra)

# Compare the aperiodic exponent results across conditions
for ind, fg in enumerate(fgs):
    print("Aperiodic exponent for condition {} is {:1.4f}".format(
        ind, np.mean(fg.get_params('aperiodic_params', 'exponent'))))





files = glob.glob(
    '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/SNR/AfterWhole/ICAwholeClean_homogenize/1*_2*_*_Rem_rerefwhole_ICA_icapru.set')
outpath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/ASSR_fooof/'
imageoutpath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/ASSR_fooof/images/'

for fname in files[211:299]:
    sub = fname[111:125]
    savefile = outpath + sub + '_ExpOff_40Hz.csv'
    if path.exists(savefile):
        print(f"skipping {sub}")
        continue

    raw = mne.io.read_raw_eeglab(fname, preload=True)
    raw = raw.pick_types(meg=False, eeg=True, eog=False, exclude='bads')

    #raw.plot()
    raw._data = check_nans(raw._data, nan_policy='zero')
    boundary_idx = np.where(raw.annotations.description == 'boundary')[
        0]  # find where the code has triggers labeled boundary

    raw.annotations.delete(boundary_idx)  # delete boundary triggers

    raw.notch_filter(np.arange(58, 62, 2), filter_length='auto',
                     phase='zero')  # filter out the 60Hz artifact from power line noise
    # raw.plot_psd(area_mode='range', tmax=10.0, average=False);

    events_from_annot, event_dict = mne.events_from_annotations(raw)

    epochs = mne.Epochs(raw, events_from_annot, tmin=0, tmax=0.5, event_id=event_dict, preload=True, baseline=(0, 0))
    if '4' in epochs.event_id:
        evoked = epochs['4']
        evokedAvg = epochs['4'].average()
        spectrum = evokedAvg.compute_psd(method='welch', fmin=0, fmax=75, n_fft=evokedAvg.times.size)  # find spectrum

        psds, myfreqs = spectrum.get_data(return_freqs=True)  # grab frequency values corresponding to spectrum
        myspectra = np.array(spectrum._data)  # grab spectra values

        # Initialize a model object for spectral parameterization, with some settings

        subDF = pd.DataFrame()
        subPeaksDF = pd.DataFrame()

        # Extract data for the current channel

        fg = SpectralGroupModel(peak_width_limits=[0.5, 7], min_peak_height=0, max_n_peaks=7, verbose=False,
                                peak_threshold=0.5)

        fg.fit(myfreqs, myspectra, [3, 55])

        for channel in range(myspectra.shape[1]):
            fm = fg.get_model(ind=channel, regenerate=True)
            # # # # # Print result and plot extracted model fit
            fm.print_results()
            fm.plot()
            plt.show()
            savefile = imageoutpath + sub + '_' + raw.ch_names[channel] + '_40hz'
            plt.savefig(savefile)
            plt.close()

        bands = Bands({'theta': [4, 8],
                       'alpha': [8, 13],
                       'beta': [13, 30],
                       'gamma': [30, 55]})

        # gammas = get_band_peak(fm, bands.gamma)

        # Extract aperiodic parameters
        Aps = fg.get_params('aperiodic_params')  # offset, exponent

        # Extract peak parameters
        channelPeaks = fg.get_params(
            'peak_params')  # Center freq, Power of peak above aperiodic component, bandwidth of peak
        CF = channelPeaks[:, 0]
        peakpower = channelPeaks[:, 1]
        bandwidth = channelPeaks[:, 2]

        # Extract goodness-of-fit metrics
        Errors = fg.get_params('error')
        R2s = fg.get_params('r_squared')

        subDF = pd.DataFrame(
            {'Subject': sub, 'Channel': raw.ch_names, 'Offset': Aps[:, 0], 'Exponent': Aps[:, 1], 'Error': Errors,
             'Rsquared': R2s})

        saveFile = outpath + sub + '_ExpOff_40Hz.csv'
        subDF.to_csv(saveFile)

        subPeaksDF = pd.DataFrame({'Subject': sub, 'CenterFreq': CF, 'Power': peakpower, 'Bandwidth': bandwidth,
                                   'ChannelPeaks': channelPeaks[:, 3]})

        saveFile = outpath + sub + '_Peaks_40Hz.csv'
        subPeaksDF.to_csv(saveFile)

    else:

        continue


## run fooof on the entire thing of 40hz, including rest blocks

files = glob.glob(
    '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/SNR/AfterWhole/ICAwholeClean_homogenize/1*_2*_*_Rem_rerefwhole_ICA_icapru.set')
outpath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/ASSR_fooof/'
imageoutpath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/ASSR_fooof/images_continuousData/'

for fname in files[156:299]:
    sub = fname[111:125]
    savefile = outpath + sub + '_ExpOff_40Hz.csv'
    if path.exists(savefile):
        print(f"skipping {sub}")
        continue

    raw = mne.io.read_raw_eeglab(fname, preload=True)
    raw = raw.pick_types(meg=False, eeg=True, eog=False, exclude='bads')

    #raw.plot()
    raw._data = check_nans(raw._data, nan_policy='zero')
    boundary_idx = np.where(raw.annotations.description == 'boundary')[
        0]  # find where the code has triggers labeled boundary

    raw.annotations.delete(boundary_idx)  # delete boundary triggers

    raw.notch_filter(np.arange(58, 62, 2), filter_length='auto',
                     phase='zero')  # filter out the 60Hz artifact from power line noise
    # raw.plot_psd(area_mode='range', tmax=10.0, average=False);

    events_from_annot, event_dict = mne.events_from_annotations(raw)

    events_to_remove = events_from_annot[np.isin(events_from_annot[:, -1], [2, 3])]

    stim40_idx = np.where((raw.annotations.description != '4'))[0]  # find indices where the trigger DOES NOT equal 16129, 15361, 0(eyes closed) or 16130, 15262, 1 (eyes open)
    closedData = raw.copy()  # copy the raw file
    closedData.annotations.delete(stim40_idx)  # delete all other triggers

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

    spectrum = appendedRawSegments_editBoundary.compute_psd(method='welch', fmin=0, fmax=75, n_fft=512)  # find spectrum

    psds, myfreqs = spectrum.get_data(return_freqs=True)  # grab frequency values corresponding to spectrum
    myspectra = np.array(spectrum._data)  # grab spectra values

    # Initialize a model object for spectral parameterization, with some settings

    subDF = pd.DataFrame()
    subPeaksDF = pd.DataFrame()

    # Extract data for the current channel

    fg = SpectralGroupModel(min_peak_height=0, max_n_peaks=7, verbose=False)

    fg.fit(myfreqs, myspectra, [3, 55])
    fm = fg.get_model(ind=23, regenerate=True)
    # # # # # Print result and plot extracted model fit
    fm.print_results()
    fm.plot()
    plt.show()

    for channel in range(myspectra.shape[0]):
        fm = fg.get_model(ind=channel, regenerate=True)
        # # # # # Print result and plot extracted model fit
        fm.print_results()
        fm.plot()
        plt.show()
        savefile = imageoutpath + sub + '_' + raw.ch_names[channel] + '_40hz'
        plt.savefig(savefile)
        plt.close()

    bands = Bands({'theta': [4, 8],
                   'alpha': [8, 13],
                   'beta': [13, 30],
                   'gamma': [30, 55]})

    #gammas = get_band_peak(fm, bands.gamma)

    # Extract aperiodic parameters
    Aps = fg.get_params('aperiodic_params')  # offset, exponent

    # Extract peak parameters
    channelPeaks = fg.get_params(
        'peak_params')  # Center freq, Power of peak above aperiodic component, bandwidth of peak
    CF = channelPeaks[:, 0]
    peakpower = channelPeaks[:, 1]
    bandwidth = channelPeaks[:, 2]

    # Extract goodness-of-fit metrics
    Errors = fg.get_params('error')
    R2s = fg.get_params('r_squared')

    subDF = pd.DataFrame(
        {'Subject': sub, 'Channel': raw.ch_names, 'Offset': Aps[:, 0], 'Exponent': Aps[:, 1], 'Error': Errors,
         'Rsquared': R2s})

    saveFile = outpath + sub + '_ExpOff_40Hz.csv'
    subDF.to_csv(saveFile)

    subPeaksDF = pd.DataFrame({'Subject': sub, 'CenterFreq': CF, 'Power': peakpower, 'Bandwidth': bandwidth,
                               'ChannelPeaks': channelPeaks[:, 3]})

    saveFile = outpath + sub + '_Peaks_40Hz.csv'
    subPeaksDF.to_csv(saveFile)
