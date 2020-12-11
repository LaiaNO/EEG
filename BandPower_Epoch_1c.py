import matplotlib
import matplotlib.pyplot as plt # plotting
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import scipy
import scipy.interpolate
import mne
import scipy
import statistics
from pywt import wavedec
from scipy import signal
import matplotlib.colors as mcolors
import random
import seaborn as sns
import os
from scipy import stats
import statsmodels.api as sm
import numpy as np
import matplotlib.pyplot as plt
from numpy import asarray
from numpy import savetxt

from Definicions import upperchanel
from Definicions import group_inf
from Definicions import mediumchanels
from Definicions import opteciogrups
from Definicions import find_nearest
from Definicions import grabt
from Definicions import epoch_return
from Definicions import calcul

def bandpower(data, sf, band, window_sec=None, relative=False):
    """Compute the average power of the signal x in a specific frequency band.

    Parameters
    ----------
    data : 1d-array
        Input signal in the time-domain.
    sf : float
        Sampling frequency of the data.
    band : list
        Lower and upper frequencies of the band of interest.
    window_sec : float
        Length of each window in seconds.
        If None, window_sec = (1 / min(band)) * 2
    relative : boolean
        If True, return the relative power (= divided by the total power of the signal).
        If False (default), return the absolute power.

    Return
    ------
    bp : float
        Absolute or relative band power.
    """
    from scipy.signal import welch
    from scipy.integrate import simps
    band = np.asarray(band)
    low, high = band


    # Define window length
    if window_sec is not None:
        nperseg = window_sec * sf
    else:
        nperseg = (2 / low) * sf

    # Compute the modified periodogram (Welch)
    freqs, psd = welch(data, sf, nperseg=nperseg)
    #print('freqs', freqs)
    #print('psd', psd)

    # Frequency resolution
    freq_res = freqs[1] - freqs[0]
    #print('frq',freq_res)
    
    #ok
    # Find closest indices of band in frequency vector
    idx_band = np.logical_and(freqs >= low, freqs <= high)
    #print('ind', idx_band)

    # Integral approximation of the spectrum using Simpson's rule.
    bp = simps(psd[idx_band], dx=freq_res)
    #print('bp', bp)

    if relative:
        bp /= simps(psd, dx=freq_res)
    return bp



def BanPoer_Epoch(EO_EC_Pacients, numchanel, eovsEO):

    pacient_beta_EO = []
    pacient_alpha_EO = []
    pacient_gama_EO = []
    pacient_theta_EO = []
    pacient_delta_EO = []

    chanels_betta = []
    chanels_gama = []
    chanels_alpha = []
    chanels_theta = []
    chanels_delta = []

    for pacient in EO_EC_Pacients:
        EO_Epochs = epoch_return(pacient[eovsEO])

        for epoch in EO_Epochs[numchanel]:
            f_beta = bandpower(epoch, 250,[30, 250], window_sec=40, relative=True)
            chanels_betta.append(f_beta)
            #gamma
            f_gama = bandpower(epoch, 250, [32, 100], window_sec=40, relative=True)
            chanels_gama.append(f_gama)
            #alpha
            f_alpha = bandpower(epoch, 250, [9, 13], window_sec=40, relative=True)
            chanels_alpha.append(f_alpha)
            #theta
            f_theta = bandpower(epoch, 250, [4, 8], window_sec=40, relative=True)
            chanels_theta.append(f_theta)
            #delta
            f_delta = bandpower(epoch, 250, [0.1, 4], window_sec=40, relative=True)
            chanels_delta.append(f_delta)
                

    mean_x = statistics.median(chanels_betta)
    pacient_beta_EO.append(mean_x)

    mean_x = statistics.median(chanels_gama)
    pacient_alpha_EO.append(mean_x)

    mean_x = statistics.median(chanels_alpha)
    pacient_gama_EO.append(mean_x)

    mean_x = statistics.median(chanels_theta)
    pacient_theta_EO.append(mean_x)

    mean_x = statistics.median(chanels_delta)
    pacient_delta_EO.append(mean_x)
    return pacient_beta_EO, pacient_alpha_EO, pacient_gama_EO, pacient_theta_EO, pacient_delta_EO
