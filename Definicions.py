
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

'''FOR THE INITIAL DATA TREATMENT'''
#DEF UPPER CHANNEL
def upperchanel(chanlesnames):
    upperchanels=[]
    for i in chanlesnames:
        r = i.upper()
        upperchanels.append(r)
    return upperchanels

#DEF GRUP INFO DE CADA COMPONENT
def group_inf(nomgroup, data_chan, upperchanels):
    group_date = []

    for e in nomgroup:
        if e in upperchanels:
            num = upperchanels.index(e)
            #print(num)
            group_date.append(data_chan[num])  
    #print(group_date)
    return group_date

#DEF MEDIUM DELS GRUPS
#group_dat = [[1,2,3,..., 119345],[1,2,3,..., 119345],[1,2,3,..., 119345],[1,2,3,..., 119345]]
def mediumchanels(group_dat):
    segundo = len(group_dat[0])
    primero = len(group_dat)
    mean_data = []
    for i in range (0,segundo):
        suma1=[]
        for e in range(0,primero):
            suma1.append(group_dat[e][i])
        numpero = statistics.median(suma1)
        mean_data.append(numpero)
    return mean_data
def opteciogrups(chanlesnames, groupp, data_chan):
    group_date_finale = [] #Save all informatio chanels 12.
    group_date_final = [] #final copy
    upchan = upperchanel(chanlesnames)
    for i in groupp:
        group_date = group_inf(i, data_chan, upchan)
        mean_group = mediumchanels(group_date)
        mean_array = np.array(mean_group)
        group_date_finale.append(mean_array)
    group_date_final = np.array(group_date_finale)
    return group_date_final



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
    
def epoch_return(groups_date_finalle2):
    epoch_dat=[]
    for data in groups_date_finalle2:
        Fs = 250.0
        Ts = 1.0/Fs
        t = np.arange(len(data)) / Fs

        n = len(data) # length of the signal
        k = np.arange(n)
        T = n/Fs
        frq = k/T # two sides frequency range
        frq = frq[range(int(n/2))]

        Y = np.fft.fft(data)/n
        Y = Y[range(int(n/2))]
        save_epoc_data = []
        save_epoc_temps = []
        tinici = 0
        tfinal = 1250
        for i in range(0,int(int(T)/5)):
            #FROM THE CHANEL SELECTED EXTRACT THE t DATA
            t = grabt(data)

            #SELECT FROM t AND group_date THE DATA FROM POSITIONS tinici to tfinal
            epoc_temps = t[tinici:tfinal]
            epoc_data = data[tinici:tfinal]

            save_epoc_data.append(epoc_data)
            save_epoc_temps.append(epoc_temps)    

            tfinal=tfinal+1250
            tinici= tinici+1250
        epoch_dat.append(save_epoc_data)
    return epoch_dat

'''def collect(l, index):
   return map(itemgetter(index), l)'''

def calcul(names, sorted_list_EO, TeiQueSF_emotionality, pacient_beta, yo, oy, color):
    all_info = []
    x_corr = []
    y_corr = []
    for i in range(0, len(sorted_list_EO)):
        hename = sorted_list_EO[i]
        hename = str(hename[:-7])
        if hename in names:
            indices = [i for i, s in enumerate(names) if hename in s]
            x = (float(TeiQueSF_emotionality[int(indices[0])]))
            x_corr.append(x)
            y = (float(pacient_beta[i]))
            y_corr.append(y)
            plt.scatter(x,y, c=color)
            
    #LINE EC
    z = np.poly1d(np.polyfit(x_corr, y_corr, 1))
    y_len = np.array(len(y_corr))
    xp = np.linspace(yo, oy, y_len)
    y = z(xp)
    plt.plot(xp, y, c='b')

    #PLOT ALL POINTS
    plt.title('Plot show the correlation between TeiQueSF_emotionality and Band Power')
    print(stats.pearsonr(x_corr, y_corr))

    return stats.pearsonr(x_corr, y_corr)