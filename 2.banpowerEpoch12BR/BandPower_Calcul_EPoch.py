import matplotlib
import matplotlib.pyplot as plt  # plotting
import numpy as np  # linear algebra
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
import scipy
import scipy.interpolate
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
import scipy.cluster

from scipy.signal import welch
from scipy.integrate import simps

'''FOR THE EPOCH definiton.'''


def epoch_return(chanels, numchanel):
    """ Return the data devided by epochs of 5 sec.
    Parameters
    ---------- 
    chanels = the 12 brain reagions; Len chanels[0/11] = 120000 - 104000 (domain amplitude)
    numchanel = channel number [0/11]
    # Other features not included in the analisis; Ts = 1.0/Fs  # Sampling Time // t = np.arange(n)/Fs  # Time  // #T = t/n
    """
    data = chanels[numchanel]  # Select the data of the chanel
    Fs = 250.0  # Our frequency sample
    n = len(data)  # length of the signal
    T = n/Fs  # Period in time domain

    save_epoc_data = []
    tinici = 0
    tfinal = 1250
    for i in range(0, int(int(T)/5)):
        # SELECT FROM group_date THE DATA FROM POSITIONS tinici to tfinal
        epoc_data = data[tinici:tfinal]
        save_epoc_data.append(epoc_data)

        tfinal = tfinal+1250
        tinici = tinici+1250
    return save_epoc_data


'''BANDPOWER EXTRACTION'''


def bandpower(data, sf, band, window_sec=None, relative=False):
    """ SOURCE: https://raphaelvallat.com/bandpower.html
    Compute the average power of the signal x in a specific frequency band.
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

    band = np.asarray(band)
    low, high = band

    # Define window length
    if window_sec is not None:
        nperseg = window_sec * sf
    else:
        nperseg = (2 / low) * sf

    # Compute the modified periodogram (Welch)
    freqs, psd = welch(data, sf, nperseg=nperseg)

    freq_res = freqs[1] - freqs[0]  # Frequency resolution

    # Find closest indices of band in frequency vector
    idx_band = np.logical_and(freqs >= low, freqs <= high)

    # Integral approximation of the spectrum using Simpson's rule.
    bp = simps(psd[idx_band], dx=freq_res)

    if relative:
        bp /= simps(psd, dx=freq_res)
    return bp


'''Execute all'''


def BanPoer_Epoch(EO_EC_Data, numchanel, EOorEC):
    """ Return the BandPower for each epoch for each Band.
    Parameters
    ---------- 
    EO_EC_Data:
        List of all the patients, each one with EO and EC, with their num. of chanels reduced to 12. Ordered by name
        Len = 187 (num. of patients)
        Len EO_EC_Data[x] = 2
        Len EO_EC_Data[x][0/1] = 12
        Len EO_EC_Data[x][0/1][0/11] = 120000 - 104000
    numchanel: The Brain Ragions selected of each patient. 
    EOorEC: Select if we what the EO or the EC data from EO_EC_Data
    """
    pacient_beta_EO = []
    pacient_alpha_EO = []
    pacient_gama_l_EO = []
    pacient_theta_EO = []
    pacient_delta_EO = []

    returnEpochBR = []  # SAVE ALL THE BANDS

    for pacient in EO_EC_Data:
        EO_Epochs = epoch_return(pacient[EOorEC], numchanel)
        chanels_betta = []
        chanels_gama_l = []
        chanels_alpha = []
        chanels_theta = []
        chanels_delta = []

        for epoch in EO_Epochs:
            # gamma lower
            f_gama_lower = bandpower(
                epoch, 250, [30, 45], window_sec=1, relative=True)
            chanels_gama_l.append(f_gama_lower)
            # beta
            f_beta = bandpower(
                epoch, 250, [12, 30], window_sec=1, relative=True)
            chanels_betta.append(f_beta)
            # alpha
            f_alpha = bandpower(
                epoch, 250, [8, 12], window_sec=1, relative=True)
            chanels_alpha.append(f_alpha)
            # theta
            f_theta = bandpower(
                epoch, 250, [4, 8], window_sec=1, relative=True)
            chanels_theta.append(f_theta)
            # delta
            f_delta = bandpower(
                epoch, 250, [2, 4], window_sec=1, relative=True)
            chanels_delta.append(f_delta)

        mean_x = statistics.median(chanels_betta)
        mean_l = statistics.median(chanels_gama_l)
        mean_t = statistics.median(chanels_alpha)
        mean_o = statistics.median(chanels_theta)
        mean_p = statistics.median(chanels_delta)

        pacient_beta_EO.append(mean_x)
        pacient_alpha_EO.append(mean_t)
        pacient_gama_l_EO.append(mean_l)
        pacient_theta_EO.append(mean_o)
        pacient_delta_EO.append(mean_p)

    returnEpochBR.append(pacient_delta_EO)
    returnEpochBR.append(pacient_theta_EO)
    returnEpochBR.append(pacient_alpha_EO)
    returnEpochBR.append(pacient_beta_EO)
    returnEpochBR.append(pacient_gama_l_EO)
    return returnEpochBR
