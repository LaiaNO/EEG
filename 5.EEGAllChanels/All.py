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
from numpy import asarray
from numpy import savetxt
import pickle
import seaborn as sns
import os
from scipy import stats
import statsmodels.api as sm
import scipy.cluster
import scipy
import numpy as np
import matplotlib.pyplot as plt
from numpy import asarray
from numpy import savetxt

def Import_Patients(basePATH):

    ###READ FILES IN THE FOLDER AND SAVE THEM

    #Save the list of names separated EO vs EC
    listPATHEC = []
    listPATHEO = []
    #Save only the set files given that in the upload they automaticly search for the .fdt
    for root, dirs, files in os.walk(str(basePATH)):
        for filename in files:
            if '.set' in filename:
                if 'EC' in filename:
                    listPATHEC.append(filename)
                if 'EO' in filename:
                    listPATHEO.append(filename)
    #Save again but sorted by the name that way we know they are in order
    sorted_list_EO = sorted(listPATHEO)
    sorted_list_EC = sorted(listPATHEC)



    #ONCE WE HAVE THE NAMES IS TIME TO UPLOAD THE DATA
    #Variable were we will have all the patients with their EO and EC respective.
    EO_EC_Pacients = []
    Chanels = []

    for i in sorted_list_EO:
        #Load the doc
        x=mne.io.read_raw_eeglab(basePATH+i, preload=True, verbose=True)
        #GET DATA
        EO = x._data
        #GET CHANELS
        chanles_names_EO = x.ch_names

        #Select the same but with EC
        EC_name = [sorted_list_EC.index(e) for e in sorted_list_EC if i[:-6] in e]
        EC_name_PATH = sorted_list_EC[EC_name[0]]
        #LOAD EC
        x2=mne.io.read_raw_eeglab(basePATH+EC_name_PATH, preload=True, verbose=True)
        #GET DATA
        EC = x2._data
        #GET CHANELS
        chanles_names_EC = x2.ch_names

        #save
        EO_EC_P = []
        EO_EC_P.append(EO)
        EO_EC_P.append(EC)
        EO_EC_Pacients.append(EO_EC_P)
        chan = []
        chan.append(chanles_names_EO)
        chan.append(chanles_names_EC)
        Chanels.append(chan)
    
    return(sorted_list_EO, sorted_list_EC, EO_EC_Pacients, Chanels)


#band
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

# x=datos
# fs=sample frequency
# fmin=min frequency
# fmax=max frequency


def BanPoer_Epoch(EO_EC_Pacients, eovsEO):

    beta_EO = []
    alpha_EO = []
    gama_l_EO = []
    gama_u_EO = []
    theta_EO = []
    delta_EO = []


    for pacient in EO_EC_Pacients:
        pacient_beta_EO = []
        pacient_alpha_EO = []
        pacient_gama_l_EO = []
        pacient_gama_u_EO = []
        pacient_theta_EO = []
        pacient_delta_EO = []
        for e in range(0,len(pacient[eovsEO])):
            EO_Epochs = epoch_return(pacient[eovsEO], e)

            chanels_betta = []
            chanels_gama_u = []
            chanels_gama_l = []
            chanels_alpha = []
            chanels_theta = []
            chanels_delta = []

            for epoch in EO_Epochs:
            
                #gamma upper
                f_gama_upper = bandpower(epoch, 250, [45, 250], window_sec=1, relative=True)
                chanels_gama_u.append(f_gama_upper)
                #gamma lower
                f_gama_lower = bandpower(epoch, 250, [30, 45], window_sec=1, relative=True)
                chanels_gama_l.append(f_gama_lower)
                #beta
                f_beta = bandpower(epoch, 250,[12, 30], window_sec=1, relative=True)
                chanels_betta.append(f_beta)
                #alpha
                f_alpha = bandpower(epoch, 250, [8, 12], window_sec=1, relative=True)
                chanels_alpha.append(f_alpha)
                #theta
                f_theta = bandpower(epoch, 250, [4, 8], window_sec=1, relative=True)
                chanels_theta.append(f_theta)
                #delta
                f_delta = bandpower(epoch, 250, [2, 4], window_sec=1, relative=True)
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

        beta_EO.append(pacient_beta_EO)
        alpha_EO.append(pacient_alpha_EO)
        gama_l_EO.append(pacient_gama_l_EO)
        gama_u_EO.append(pacient_gama_u_EO)
        theta_EO.append(pacient_theta_EO)
        delta_EO.append(pacient_delta_EO)

    return delta_EO, theta_EO, alpha_EO, beta_EO, gama_l_EO 

'''FOR THE EPOCH'''
#TO FIND THE NEAREST POSITION OF THE DATA THAT CORESPOND TO THE ONE INTRODUCED.
def find_nearest(array,value): 
       idx = (np.abs(array-value)).argmin()
       return int(idx)

#EXTRACT THE t FROM THE DATA
def grabt(dat):
    Fs = 250.0
    Ts = 1.0/Fs
    t = np.arange(len(dat)) / Fs
    return t

#  groups_date_finalle2=canales
#  numchanel=channel number
def epoch_return(groups_date_finalle2, numchanel):
    epoch_dat=[]
    data = groups_date_finalle2[numchanel]
    Fs = 250.0
    # Ts = 1.0/Fs
    n = len(data) # length of the signal
    t = np.arange(n) / Fs

    k = np.arange(n)
    T = n/Fs
    #T= t/n

    save_epoc_data = []
    # save_epoc_temps = []
    tinici = 0
    tfinal = 1250
    for i in range(0,int(int(T)/5)):
    #for i in range(0,int(int(t)/5)):
        #FROM THE CHANEL SELECTED EXTRACT THE t DATA
        #t = grabt(data)

        #SELECT FROM t AND group_date THE DATA FROM POSITIONS tinici to tfinal
        #epoc_temps = t[tinici:tfinal]
        epoc_data = data[tinici:tfinal]

        save_epoc_data.append(epoc_data)
        #save_epoc_temps.append(epoc_temps)    

        tfinal=tfinal+1250
        tinici= tinici+1250
    #epoch_dat.append(save_epoc_data)
    return save_epoc_data

'''def collect(l, index):
   return map(itemgetter(index), l)'''



def epoch_return_todo(groups_date_finalle2, numchanel):
    epoch_dat=[]
    for data in groups_date_finalle2:
        Fs = 250.0
        Ts = 1.0/Fs
        t = np.arange(len(data)) / Fs

        n = len(data) # length of the signal
        k = np.arange(n)
        T = n/Fs
        #frq = k/T # two sides frequency range
        #frq = frq[range(int(n/2))]

        #Y = np.fft.fft(data)/n
        #Y = Y[range(int(n/2))]

        save_epoc_data = []
        save_epoc_temps = []
        tinici = 0
        tfinal = 1250
        for i in range(0,int(int(T)/5)):
            #FROM THE CHANEL SELECTED EXTRACT THE t DATA
            t = grabt(data)

            #SELECT FROM t AND group_date THE DATA FROM POSITIONS tinici to tfinal
            #epoc_temps = t[tinici:tfinal]
            epoc_data = data[tinici:tfinal]

            save_epoc_data.append(epoc_data)
            #save_epoc_temps.append(epoc_temps)    

            tfinal=tfinal+1250
            tinici= tinici+1250
        #epoch_dat.append(save_epoc_data)
    return save_epoc_data

'''def collect(l, index):
   return map(itemgetter(index), l)'''


#UPLOAD ALL THE FILES OF THE PATIENS IN A LIST LIKE THE FOLLOWING:
#Subject = EO_EC_Pacients[S] 
#EO = EO_EC_Pacients[S][0]
#EC = EO_EC_Pacients[S][1]
#Channel = EO_EC_Pacients[0][1][Num] 

#Subject = Chanels[S] 
#EO_Chanels_Names = Chanels[S][0]
#EC_Chanels_Names = Chanels[S][1]
#Channel_Name = Chanels[0][1][Num] 
basePATH = '/Users/laianavarroolivella/Proyectos/EEG/Files_EEG/Preprocessed/'
sorted_list_EO, sorted_list_EC, EO_EC_Pacients, Chanels = Import_Patients(basePATH)
print(len(sorted_list_EO))


p1O, p2O, p3O, p4O, p5O = BanPoer_Epoch(EO_EC_Pacients, 0)
p1C, p2C, p3C, p4C, p5C = BanPoer_Epoch(EO_EC_Pacients, 1)
print(len(p1C))
listageneralO = []
listageneralC = []

listageneralO.append(p2O)
listageneralO.append(p3O)
listageneralO.append(p4O)
listageneralO.append(p5O)
print(len(listageneralO))
listageneralC.append(p2C)
listageneralC.append(p3C)
listageneralC.append(p4C)
listageneralC.append(p5C)

list_finalC = []
list_finalO = []
list_finalC.append(listageneralC)
list_finalO.append(listageneralO)

with open("/Users/laianavarroolivella/Proyectos/EEG/list_finalC.txt", "wb") as fp:   #Pickling
    pickle.dump(list_finalC, fp)

with open("/Users/laianavarroolivella/Proyectos/EEG/list_finalO.txt", "wb") as fp:   #Pickling
    pickle.dump(list_finalO, fp)


with open("/Users/laianavarroolivella/Proyectos/EEG/Chanels.txt", "wb") as fp:   #Pickling
    pickle.dump(Chanels, fp)

print("Finish")