
import matplotlib
import matplotlib.pyplot as plt # plotting
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import scipy
import scipy.interpolate
#import mne
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

def calcul(names, sorted_list_EO, TeiQueSF_emotionality, age, pacient_beta, yo, oy):
    savecor_1r = []
    savecor_2r = []
    savecor_3r = []
    
    for e in range(0,len(pacient_beta)):
        x_corr_1r = []
        y_corr_1r = []
        x_corr_2r = []
        y_corr_2r = []
        x_corr_3r = []
        y_corr_3r = []
        for i in range(0, len(sorted_list_EO)):
            hename = sorted_list_EO[i]
            hename = str(hename[:-7])
            if hename in names:
                indices = [i for i, s in enumerate(names) if hename in s]
                x = (float(TeiQueSF_emotionality[int(indices[0])]))
                g = ((age[(indices[0])]))
                if str(g)=='20-25' or str(g)=='25-30' or str(g)=='30-35' or str(g)=='35-40':
                    x_corr_1r.append(x)
                    y = (float(pacient_beta[e][i]))
                    y_corr_1r.append(y)
                if str(g)=='55-60' or str(g)=='60-65':
                    x_corr_2r.append(x)
                    y = (float(pacient_beta[e][i]))
                    y_corr_2r.append(y)
                if str(g)=='65-70' or str(g)=='70-75' or str(g)=='75-80':
                    x_corr_3r.append(x)
                    y = (float(pacient_beta[e][i]))
                    y_corr_3r.append(y)
                
        savecor_1r.append(stats.pearsonr(x_corr_1r, y_corr_1r))
        savecor_2r.append(stats.pearsonr(x_corr_2r, y_corr_2r))
        savecor_3r.append(stats.pearsonr(x_corr_3r, y_corr_3r))

    return savecor_1r, savecor_2r, savecor_3r


def minus_return_bandpower1(names, sorted_list_EO, TeiQueSF_emotionality, gender, list_bands):
    minus_5_1 = []
    minus_5_2 = []
    minus_5_3 = []
    for e in list_bands:
        y_corr_1 = []
        y_corr_2 = []
        y_corr_3 = []
        for i in range(0, len(sorted_list_EO)):
            hename = sorted_list_EO[i]
            hename = str(hename[:-7])
            if hename in names:
                indices = [i for i, s in enumerate(names) if hename in s]
                x = (float(TeiQueSF_emotionality[int(indices[0])]))
                g = ((gender[(indices[0])]))
                if str(g)=='20-25' or str(g)=='25-30' or str(g)=='30-35' or str(g)=='35-40' and x<4.0:
                    y = (float(e[i]))
                    y_corr_1.append(y)
                if str(g)=='55-60' or str(g)=='60-65' and x<4.0:
                    y = (float(e[i]))
                    y_corr_2.append(y)
                if str(g)=='65-70' or str(g)=='70-75' or str(g)=='75-80' and x<4.0:
                    y = (float(e[i]))
                    y_corr_3.append(y)

        #PLOT ALL POINTS
        if len(y_corr_1)>0:
            resultat_1 = statistics.median(y_corr_1)
            minus_5_1.append(resultat_1)
        if len(y_corr_2)>0:
            resultat_2 = statistics.median(y_corr_2)
            minus_5_2.append(resultat_2)
        if len(y_corr_3)>0:
            resultat_3 = statistics.median(y_corr_3)
            minus_5_3.append(resultat_3)
        
    return minus_5_1, minus_5_2, minus_5_3




def max5_return_bandpower1(names, sorted_list_EO, TeiQueSF_emotionality, gender, list_bands):
    minus_5_1 = []
    minus_5_2 = []
    minus_5_3 = []
    for e in list_bands:
        y_corr_1 = []
        y_corr_2 = []
        y_corr_3 = []
        for i in range(0, len(sorted_list_EO)):
            hename = sorted_list_EO[i]
            hename = str(hename[:-7])
            if hename in names:
                indices = [i for i, s in enumerate(names) if hename in s]
                x = (float(TeiQueSF_emotionality[int(indices[0])]))
                g = ((gender[(indices[0])]))
                if str(g)=='20-25' or str(g)=='25-30' or str(g)=='30-35' or str(g)=='35-40' and x>4.0:
                    y = (float(e[i]))
                    y_corr_1.append(y)
                if str(g)=='55-60' or str(g)=='60-65' and x>4.0:
                    y = (float(e[i]))
                    y_corr_2.append(y)
                if str(g)=='65-70' or str(g)=='70-75' or str(g)=='75-80' and x>4.0:
                    y = (float(e[i]))
                    y_corr_3.append(y)

        #PLOT ALL POINTS
        if len(y_corr_1)>0:
            resultat_1 = statistics.median(y_corr_1)
            minus_5_1.append(resultat_1)
        if len(y_corr_1)>0:
            resultat_2 = statistics.median(y_corr_2)
            minus_5_2.append(resultat_2)
        if len(y_corr_1)>0:
            resultat_3 = statistics.median(y_corr_3)
            minus_5_3.append(resultat_3)
        
    return minus_5_1, minus_5_2, minus_5_3


def plotbar_info_bandpower1 (names, sorted_list_EO, TeiQueSF_emotionality, gender, list_bands):
    #chanels1_delta_EO, chanels1_theta_EO, chanels1_alpha_EO, chanels1_beta_EO, chanels1_gama_l_EO, chanels1_gama_u_EO
    labels = ['Delta', 'Theta', 'Alpha', 'Beta', 'Gama Lower']
    minus5_1, minus5_2, minus5_3 = minus_return_bandpower1(names, sorted_list_EO, TeiQueSF_emotionality, gender,list_bands)
    max_5_1, max_5_2, max_5_3 = max5_return_bandpower1(names, sorted_list_EO, TeiQueSF_emotionality, gender, list_bands)
    #print(minus5)
    str1Max, str2Max, str3Max = max5_return_str1(names, sorted_list_EO, TeiQueSF_emotionality, gender, list_bands)
    str1Min, str2Min, str3Min = min5_return_str1(names, sorted_list_EO, TeiQueSF_emotionality, gender, list_bands)

    x = np.arange(len(labels))  # the label locations
    width = 0.10  # the width of the bars


    fig, ax = plt.subplots()
    if len(minus5_1)>0:
        rects1 = ax.bar(x-0.3, minus5_1, width, yerr=str1Min,capsize=0.05,label='-4 1')
    if len(max_5_1)>0:
        rects3 = ax.bar(x-0.2, max_5_1, width, yerr=str1Max,capsize=0.05, label='+4 1')
    if len(minus5_2)>0:
        rects4 = ax.bar(x-0.1, minus5_2, width, yerr=str2Min,capsize=0.05,label='-4 2')
    if len(max_5_2)>0:
        rects5 = ax.bar(x, max_5_2, width, yerr=str2Max,capsize=0.05, label='+4 2')
    if len(minus5_3)>0:
        rects2 = ax.bar(x+0.1, minus5_3, width, yerr=str3Min,capsize=0.05,label='-4 3')
    if len(max_5_3)>0:
        rects4 = ax.bar(x+0.2, max_5_3, width, yerr=str3Max,capsize=0.05, label='+4 3')


    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Band Power')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    fig.tight_layout()

    plt.show()


  

    return max_5_1, max_5_2, max_5_3





def max5_return_str1(names, sorted_list_EO, TeiQueSF_emotionality, age, list_bands):
    stt_1 = []
    stt_2 = []
    stt_3 = []
    for e in list_bands:
        y_corr_1r = []
        y_corr_2r = []
        y_corr_3r = []
        for i in range(0, len(sorted_list_EO)):
            hename = sorted_list_EO[i]
            hename = str(hename[:-7])
            if hename in names:
                indices = [i for i, s in enumerate(names) if hename in s]
                x = (float(TeiQueSF_emotionality[int(indices[0])]))
                g = ((age[(indices[0])]))
                if str(g)=='20-25' or str(g)=='25-30' or str(g)=='30-35' or str(g)=='35-40' and x>4.0:
                    y = (float(e[i]))
                    y_corr_1r.append(y)
                if str(g)=='55-60' or str(g)=='60-65' and x>4.0:
                    y = (float(e[i]))
                    y_corr_2r.append(y)
                if str(g)=='65-70' or str(g)=='70-75' or str(g)=='75-80' and x>4.0:
                    y = (float(e[i]))
                    y_corr_3r.append(y)
                
        #PLOT ALL POINTS
        resultat_m = scipy.stats.sem(y_corr_1r)
        stt_1.append(resultat_m)

        resultat_f = scipy.stats.sem(y_corr_2r)
        stt_2.append(resultat_f)

        resultat_f = scipy.stats.sem(y_corr_3r)
        stt_3.append(resultat_f)


    return stt_1, stt_2, stt_3

def min5_return_str1(names, sorted_list_EO, TeiQueSF_emotionality, age, list_bands):
    stt_1 = []
    stt_2 = []
    stt_3 = []
    for e in list_bands:
        y_corr_1r = []
        y_corr_2r = []
        y_corr_3r = []
        for i in range(0, len(sorted_list_EO)):
            hename = sorted_list_EO[i]
            hename = str(hename[:-7])
            if hename in names:
                indices = [i for i, s in enumerate(names) if hename in s]
                x = (float(TeiQueSF_emotionality[int(indices[0])]))
                g = ((age[(indices[0])]))
                if str(g)=='20-25' or str(g)=='25-30' or str(g)=='30-35' or str(g)=='35-40' and x<4.0:
                    y = (float(e[i]))
                    y_corr_1r.append(y)
                if str(g)=='55-60' or str(g)=='60-65' and x<4.0:
                    y = (float(e[i]))
                    y_corr_2r.append(y)
                if str(g)=='65-70' or str(g)=='70-75' or str(g)=='75-80' and x<4.0:
                    y = (float(e[i]))
                    y_corr_3r.append(y)
                
        #PLOT ALL POINTS
        resultat_m = scipy.stats.sem(y_corr_1r)
        stt_1.append(resultat_m)

        resultat_f = scipy.stats.sem(y_corr_2r)
        stt_2.append(resultat_f)

        resultat_f = scipy.stats.sem(y_corr_3r)
        stt_3.append(resultat_f)


    return stt_1, stt_2, stt_3
#standar error of the mean -> np.std(data, ddof=1) / np.sqrt(np.size(data))





def minus_return_rac(names, sorted_list_EO, TeiQueSF_emotionality, age, list_bands):
    info_chan_1 = []
    info_chan_2 = []
    info_chan_3 = []
    for e in list_bands:
        y_corr_1r = []
        y_corr_2r = []
        y_corr_3r = []
        for i in range(0, len(sorted_list_EO)):
            hename = sorted_list_EO[i]
            hename = str(hename[:-7])
            if hename in names:
                indices = [i for i, s in enumerate(names) if hename in s]
                x = (float(TeiQueSF_emotionality[int(indices[0])]))
                g = ((age[(indices[0])]))
                if str(g)=='20-25' or str(g)=='25-30' or str(g)=='30-35' or str(g)=='35-40' and x<4.0:
                    y = (float(e[i]))
                    y_corr_1r.append(y)
                if str(g)=='55-60' or str(g)=='60-65' and x<4.0:
                    y = (float(e[i]))
                    y_corr_2r.append(y)
                if str(g)=='65-70' or str(g)=='70-75' or str(g)=='75-80' and x<4.0:
                    y = (float(e[i]))
                    y_corr_3r.append(y)

        info_chan_1.append(y_corr_1r)
        info_chan_2.append(y_corr_2r)
        info_chan_3.append(y_corr_3r)
        
    return info_chan_1, info_chan_2, info_chan_3




def max5_return_rac(names, sorted_list_EO, TeiQueSF_emotionality, age, list_bands):
    info_chan_1 = []
    info_chan_2 = []
    info_chan_3 = []
    for e in list_bands:
        y_corr_1r = []
        y_corr_2r = []
        y_corr_3r = []
        for i in range(0, len(sorted_list_EO)):
            hename = sorted_list_EO[i]
            hename = str(hename[:-7])
            if hename in names:
                indices = [i for i, s in enumerate(names) if hename in s]
                x = (float(TeiQueSF_emotionality[int(indices[0])]))
                g = ((age[(indices[0])]))
                if str(g)=='20-25' or str(g)=='25-30' or str(g)=='30-35' or str(g)=='35-40' and x>4.0:
                    y = (float(e[i]))
                    y_corr_1r.append(y)
                if str(g)=='55-60' or str(g)=='60-65' and x>4.0:
                    y = (float(e[i]))
                    y_corr_2r.append(y)
                if str(g)=='65-70' or str(g)=='70-75' or str(g)=='75-80' and x>4.0:
                    y = (float(e[i]))
                    y_corr_3r.append(y)

        info_chan_1.append(y_corr_1r)
        info_chan_2.append(y_corr_2r)
        info_chan_3.append(y_corr_3r)
        
    return info_chan_1, info_chan_2, info_chan_3


def plotbar_info_rac (names, sorted_list_EO, TeiQueSF_emotionality, gender, list_bands):
    #chanels1_delta_EO, chanels1_theta_EO, chanels1_alpha_EO, chanels1_beta_EO, chanels1_gama_l_EO, chanels1_gama_u_EO
    labels = ['Delta', 'Theta', 'Alpha', 'Beta', 'Gama Lower']

    minus1, minus2, minus3 = minus_return_rac(names, sorted_list_EO, TeiQueSF_emotionality, gender, list_bands)
    max1, max2, max3 = max5_return_rac(names, sorted_list_EO, TeiQueSF_emotionality, gender, list_bands)

    resu1 = []
    for i in range(0,len(max1)):
        p = scipy.stats.ranksums(minus1[i],max1[i])
        resu1.append(p)
    print('20-40')
    print (resu1)
    
    resu2 = []
    for i in range(0,len(max2)):
        p = scipy.stats.ranksums(minus2[i],max2[i])
        resu2.append(p)
    print('40-60')
    print (resu2)

    resu3 = []
    for i in range(0,len(max3)):
        p = scipy.stats.ranksums(minus3[i], max3[i])
        resu3.append(p)
    print('40-60')
    print (resu3)
    
    return resu1, resu2, resu3