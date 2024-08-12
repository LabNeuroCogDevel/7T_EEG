import h5py
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
import fooof
from fooof import FOOOF
from fooof import FOOOFGroup
from fooof.utils import trim_spectrum, interpolate_spectrum
from fooof.plts import plot_spectra

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
import h5py

files = glob.glob(
    '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/sEEG/sEEG_rawData/P20N001/Rest/data.mat')
outpath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Results/sEEG/fooof/'

fname = files[0]

raw = scipy.io.loadmat(fname)
#raw = h5py.File(fname)
data = raw["data"]
chan_names = [x[0][0] for x in raw["info"]["chanNames"][0,0]]

info = mne.create_info(sfreq=1000, ch_names=chan_names)
simulated_raw = mne.io.RawArray(np.transpose(data), info)

#simulated_raw.notch_filter(np.arange(60, 61, 1), filter_length='auto',
                     #phase='zero', picks=chan_names)  # filter out the 60Hz artifact from power line noise
    # raw.plot_psd(area_mode='range', tmax=10.0, average=False);

spectrum = simulated_raw.compute_psd(method='welch', fmin=1, fmax=200, tmin=0, tmax=None,
                                                        picks='all', n_fft=256, n_overlap=128,
                                                        window='hamming')  # find spectrum

psds, freqs = spectrum.get_data(return_freqs=True, picks='all')  # grab frequency values corresponding to spectrum

#interp_range = [58, 62]
#freqs_int1, powers_int1 = fooof.utils.interpolate_spectrum(freqs, np.transpose(psds), interp_range)

spectra = spectrum._data  # grab spectra values

# Initialize a FOOOFGroup object, with desired settings
fg = FOOOFGroup(peak_width_limits=[0.5, 15], min_peak_height=0,
                    peak_threshold=1.5, aperiodic_mode='fixed', verbose=False)

# Define the frequency range to fit
freq_range = [1, 50]
# Fit the power spectrum model across all channels
# spectraSlice = spectra[i,:,:]
fg.fit(freqs, spectra, freq_range)
# Check the overall results of the group fits
fg.plot()
plt.show()
# # Report: fit the model, print the resulting parameters, and plot the reconstruction
nchan = data.shape[1];
nsqr = 11 #np.sqrt(nchan)//2*2
for i in range(data.shape[1]):
    ax = plt.subplot(11,11 , i+1)
    fm = fg.get_fooof(ind=118, regenerate=True)
    # # # # # Print results and plot extracted model fit
    # fm.print_results()
    fm.plot(ax=ax, add_legend=False)
    #break
plt.show()
# plt.savefig(savefile)
# plt.close()