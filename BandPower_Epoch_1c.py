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

from Definicions import find_nearest
from Definicions import grabt
from Definicions import epoch_return

def bandpower(data, sf, band, window_sec=None, relative=False):
    """Compute the average power of the signal x in a specific frequency band.

    Parameters
    ----------
    data : 1d-array
        Input signal in the time-domain.
    sf : float
        Sampling frequency of the data. (250)
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

import scipy 

def bandpower2(x, fs, fmin, fmax):
    f, Pxx = scipy.signal.periodogram(x, fs=fs)
    ind_min = scipy.argmax(f > fmin) - 1
    ind_max = scipy.argmax(f > fmax) - 1
    return scipy.trapz(Pxx[ind_min: ind_max], f[ind_min: ind_max])


def BanPoer_Epoch(EO_EC_Pacients, numchanel, eovsEO):

    pacient_beta_EO = []
    pacient_alpha_EO = []
    pacient_gama_l_EO = []
    pacient_gama_u_EO = []
    pacient_theta_EO = []
    pacient_delta_EO = []

    chanels_betta = []
    chanels_gama_u = []
    chanels_gama_l = []
    chanels_alpha = []
    chanels_theta = []
    chanels_delta = []

    for pacient in EO_EC_Pacients:
        EO_Epochs = epoch_return(pacient[eovsEO])

        for epoch in EO_Epochs[numchanel]:
            
            #gamma upper
            f_gama_upper = bandpower(epoch, 250, [80, 250], window_sec=50, relative=True)
            chanels_gama_u.append(f_gama_upper)
            #gamma lower
            f_gama_lower = bandpower(epoch, 250, [30, 80], window_sec=50, relative=True)
            chanels_gama_l.append(f_gama_lower)
            #beta
            f_beta = bandpower(epoch, 250,[15, 30], window_sec=50, relative=True)
            chanels_betta.append(f_beta)
            #alpha
            f_alpha = bandpower(epoch, 250, [8, 12], window_sec=50, relative=True)
            chanels_alpha.append(f_alpha)
            #theta
            f_theta = bandpower(epoch, 250, [4, 8], window_sec=50, relative=True)
            chanels_theta.append(f_theta)
            #delta
            f_delta = bandpower(epoch, 250, [2, 4], window_sec=50, relative=True)
            chanels_delta.append(f_delta)
                

        mean_x = statistics.median(chanels_betta)
        mean_l = statistics.median(chanels_gama_l)
        mean_y = statistics.median(chanels_gama_u)
        mean_t = statistics.median(chanels_alpha)
        mean_o = statistics.median(chanels_theta)
        mean_p = statistics.median(chanels_delta)

        pacient_beta_EO.append(mean_x)
        pacient_alpha_EO.append(mean_t)
        pacient_gama_u_EO.append(mean_y)
        pacient_gama_l_EO.append(mean_l)
        pacient_theta_EO.append(mean_o)
        pacient_delta_EO.append(mean_p)
        #chanels1_delta_EO, chanels1_theta_EO, chanels1_alpha_EO, chanels1_beta_EO, chanels1_gama_EO
    return pacient_delta_EO, pacient_theta_EO, pacient_alpha_EO, pacient_beta_EO, pacient_gama_l_EO, pacient_gama_u_EO 

def BanPoer_Epoch2(EO_EC_Pacients, numchanel, eovsEO):

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
            
            #gamma upper
            f_gama_upper = bandpower2(epoch, 250, 80, 250)
            chanels_gama_u.append(f_gama_upper)
            #gamma lower
            f_gama_lower = bandpower2(epoch, 250, 30, 80)
            chanels_gama_l.append(f_gama_lower)
            #beta
            f_beta = bandpower2(epoch, 250, 15, 30)
            chanels_betta.append(f_beta)
            #alpha
            f_alpha = bandpower2(epoch, 250, 8, 12)
            chanels_alpha.append(f_alpha)
            #theta
            f_theta = bandpower2(epoch, 250, 4, 8)
            chanels_theta.append(f_theta)
            #delta
            f_delta = bandpower2(epoch, 250, 2, 4)
            chanels_delta.append(f_delta)
                

        mean_x = statistics.median(chanels_betta)
        mean_l = statistics.median(chanels_gama_l)
        mean_y = statistics.median(chanels_gama_u)
        mean_t = statistics.median(chanels_alpha)
        mean_o = statistics.median(chanels_theta)
        mean_p = statistics.median(chanels_delta)

        pacient_beta_EO.append(mean_x)
        pacient_alpha_EO.append(mean_t)
        pacient_gama_u_EO.append(mean_y)
        pacient_gama_l_EO.append(mean_l)
        pacient_theta_EO.append(mean_o)
        pacient_delta_EO.append(mean_p)

    return pacient_delta_EO, pacient_theta_EO, pacient_alpha_EO, pacient_beta_EO, pacient_gama_l_EO, pacient_gama_u_EO 
