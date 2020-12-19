
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

#DEF




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
    
def epoch_return(groups_date_finalle2, numchanel):
    epoch_dat=[]
    data = groups_date_finalle2[numchanel]
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
