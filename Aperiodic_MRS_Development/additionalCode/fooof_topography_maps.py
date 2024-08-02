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


# Plot the topography of gamma power
chanNames = raw.ch_names
fig, axesloc = plt.subplots(nrows=1, ncols=2, figsize=(3, 4.5),
                            gridspec_kw=dict(height_ratios=[4]))

plot_topomap(gammas_pw, raw.info, cmap=cm.RdBu, contours=6, axes=axesloc[0], names=chanNames)
axesloc[1].axis('off')
divider = make_axes_locatable(axesloc[1])
cax = divider.append_axes('right', size='5%', pad=0.2)
clim = dict(kind='value', lims=[4, 8, 12])
cbar = mne.viz.plot_brain_colorbar(cax, clim, 'RdBu', label='Power')
plt.show()

# Extract alpha peaks
alphas = get_band_peak_fg(fg, bands.alpha)

# Extract the power values from the detected peaks
alpha_pw = alphas[:, 1]
# Plot the topography of alpha power
plot_topomap(alpha_pw, raw.info, cmap=cm.RdBu, contours=0);

# Plot the topographies across different frequency bands
fig, axes = plt.subplots(1, 4, figsize=(15, 5))
for ind, (label, band_def) in enumerate(bands):
    # Get the power values across channels for the current band
    band_power = check_nans(get_band_peak_fg(fg, band_def)[:, 1])

    # Create a topomap for the current oscillation band
    mne.viz.plot_topomap(band_power, raw.info, cmap=cm.viridis, contours=0,
                         axes=axes[ind], show=False);

    # Set the plot title
    axes[ind].set_title(label + ' power', {'fontsize': 20})

# plot freq by power
fig, axes = plt.subplots(1, 4, figsize=(15, 6))
for ind, (label, band_def) in enumerate(bands):
    # Get the power values across channels for the current band
    band_power = check_nans(get_band_peak_fg(fg, band_def)[:, 1])

    # Extracted and plot the power spectrum model with the most band power
    fg.get_fooof(np.argmax(band_power)).plot(ax=axes[ind], add_legend=False)

    # Set some plot aesthetics & plot title
    axes[ind].yaxis.set_ticklabels([])
    axes[ind].set_title('biggest ' + label + ' peak', {'fontsize': 16})
    plt.show()

# Extract aperiodic exponent values
exps = fg.get_params('aperiodic_params', 'exponent')
# Plot the topography of aperiodic exponents
plot_topomap(exps, raw.info, cmap=cm.viridis, contours=0)

# Compare the power spectra between low and high exponent channels
fig, ax = plt.subplots(figsize=(8, 6))
plot_spectrum(fg.freqs, fg.get_fooof(np.argmin(exps)).power_spectrum,
              ax=ax, label='Low Exponent')
plot_spectrum(fg.freqs, fg.get_fooof(np.argmax(exps)).power_spectrum,
              ax=ax, label='High Exponent')
plt.show()
