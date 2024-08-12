
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
from fooof.bands import Bands
from fooof.analysis import get_band_peak_fg
from fooof.utils import trim_spectrum
from fooof.utils.data import subsample_spectra
from fooof.sim.gen import gen_aperiodic
from fooof.data import FOOOFSettings
from fooof.plts.templates import plot_hist
from fooof.plts.spectra import plot_spectra
from fooof.plts.periodic import plot_peak_fits, plot_peak_params
from fooof.plts.aperiodic import plot_aperiodic_params, plot_aperiodic_fits

# Import functions to examine frequency-by-frequency error of model fits
from fooof.analysis.error import compute_pointwise_error_fm, compute_pointwise_error_fg

# Import helper utility to access data
from fooof.utils.download import fetch_fooof_data

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


files = glob.glob(
    '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/SNR/AfterWhole/ICAwholeClean_homogenize/1*_2*_*_Rem_rerefwhole_ICA_icapru.set')
outpath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/ASSR_fooof/'


for fname in files:
    sub = fname[111:125] + '_40Hz.npz'
    savefile = outpath + sub
    if path.exists(savefile):
        print(f"skipping {sub}")
        continue

    raw = mne.io.read_raw_eeglab(fname, preload=True)
    raw = raw.pick_types(meg=False, eeg=True, eog=False, exclude='bads')

    # raw.plot()
    raw._data = check_nans(raw._data, nan_policy='zero')
    boundary_idx = np.where(raw.annotations.description == 'boundary')[
        0]  # find where the code has triggers labeled boundary

    raw.annotations.delete(boundary_idx)  # delete boundary triggers

    raw.notch_filter(np.arange(60, 61, 1), filter_length='auto',
                     phase='zero')  # filter out the 60Hz artifact from power line noise
    # raw.plot_psd(area_mode='range', tmax=10.0, average=False);

    events_from_annot, event_dict = mne.events_from_annotations(raw)

    epochs = mne.Epochs(raw, events_from_annot, tmin=-0.2, tmax=0.5, event_id=event_dict, preload=True)

    evoked = epochs['4'].average()

    spectrum = evoked.compute_psd(method='welch', fmin=1, fmax=50, tmin=0, tmax=None,
                                                            picks='all', n_fft=256, n_overlap=128,
                                                            window='hamming')  # find spectrum
    psds, freqs = spectrum.get_data(return_freqs=True)  # grab frequency values corresponding to spectrum
    spectra = spectrum._data  # grab spectra values

    # Initialize a FOOOFGroup object, with desired settings
    fg = FOOOFGroup(peak_width_limits=[0.5,12], min_peak_height=0,
                    peak_threshold=1, aperiodic_mode='fixed', max_n_peaks=5, verbose=False)

    # Define the frequency range to fit
    freq_range = [1, 55]
    # Fit the power spectrum model across all channels
    # spectraSlice = spectra[i,:,:]
    fg.fit(freqs, spectra, freq_range)
    # Check the overall results of the group fits
    fg.plot()
    plt.show()
    # # Report: fit the model, print the resulting parameters, and plot the reconstruction
    fm = fg.get_fooof(ind=41, regenerate=True)
    # # # # # Print results and plot extracted model fit
    fm.print_results()
    fm.plot()
    # plt.ylim([-2.5, 2.5])
    #plt.show()
    # plt.savefig(savefile)
    # plt.close()

    # Extract aperiodic parameters
    aps = fg.get_params('aperiodic_params')
    exps = fg.get_params('aperiodic_params', 'exponent')

    # Extract peak parameters
    peaks = fg.get_params('peak_params')
    cfs = fg.get_params('peak_params', 'CF')

    # Extract goodness-of-fit metrics
    errors = fg.get_params('error')
    r2s = fg.get_params('r_squared')

    sub = fname[103:117] + '_40Hz'
    savefile = outpath+sub
    np.savez_compressed(savefile, aperiodic=aps, exponents=exps, peakParam=peaks, centerFreq=cfs, Error=errors, Rsquared=r2s)


bands = Bands({'theta' : [4, 8],
               'alpha' : [8, 13],
               'beta' : [13, 30],
               'gamma': [30,55]})

thetas = get_band_peak_fg(fg, bands.theta)
alphas = get_band_peak_fg(fg, bands.alpha)
betas = get_band_peak_fg(fg, bands.beta)
gammas=get_band_peak_fg(fg, bands.gamma)

plot_aperiodic_fits(aps, fg.freq_range)
plot_aperiodic_params(aps)

# Plot the parameters for peaks, split up by band
_, axes = plt.subplots(1, 4, figsize=(14, 7))
all_bands = [thetas, alphas, betas, gammas]
for ax, label, peaks in zip(np.ravel(axes), bands.labels, all_bands):
    plot_peak_params(peaks, ax=ax)
    ax.set_title(label + ' peaks', fontsize=24)
plt.subplots_adjust(hspace=0.4)

plot_peak_fits(gammas)
