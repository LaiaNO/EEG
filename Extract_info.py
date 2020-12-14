
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



def minus_return_bandpower(names, sorted_list_EO, TeiQueSF_emotionality, list_bands):
    all_info = []
    x_corr = []
    y_corr = []
    minus_5 = []
    for e in list_bands:
        for i in range(0, len(sorted_list_EO)):
            hename = sorted_list_EO[i]
            hename = str(hename[:-7])
            if hename in names:
                indices = [i for i, s in enumerate(names) if hename in s]
                x = (float(TeiQueSF_emotionality[int(indices[0])]))
                if x<5.0:
                    x_corr.append(x)
                    y = (float(e[i]))
                    y_corr.append(y)

        #PLOT ALL POINTS
        '''resultat = stats.pearsonr(x_corr, y_corr)
        minus_5.append(resultat[0])
        '''
    return y_corr




def max5_return_bandpower(names, sorted_list_EO, TeiQueSF_emotionality, list_bands):
    all_info = []
    x_corr = []
    y_corr = []
    max_5 = []
    for e in list_bands:
        for i in range(0, len(sorted_list_EO)):
            hename = sorted_list_EO[i]
            hename = str(hename[:-7])
            if hename in names:
                indices = [i for i, s in enumerate(names) if hename in s]
                x = (float(TeiQueSF_emotionality[int(indices[0])]))
                if x>5.0:
                    x_corr.append(x)
                    y = (float(e[i]))
                    y_corr.append(y)

        #PLOT ALL POINTS
        '''resultat = stats.pearsonr(x_corr, y_corr)
        max_5.append(resultat[0])'''


    return y_corr


def plotbar_info_bandpower (names, sorted_list_EO, TeiQueSF_emotionality, list_bands):
    labels = ['Betta', 'Gamma', 'ALpha', 'theta', 'delta']
    minus5 = minus_return_bandpower(names, sorted_list_EO, TeiQueSF_emotionality, list_bands)
    max5 = max5_return_bandpower(names, sorted_list_EO, TeiQueSF_emotionality, list_bands)

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, minus5, width, label='-5')
    rects2 = ax.bar(x + width/2, max5, width, label='+5')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Band Power')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()


    '''def autolabel(rects):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')

    autolabel(rects1)
    autolabel(rects2)'''

    fig.tight_layout()

    plt.show()


def minus_return_pval(names, sorted_list_EO, TeiQueSF_emotionality, list_bands):
    all_info = []
    x_corr = []
    y_corr = []
    minus_5 = []
    for e in list_bands:
        for i in range(0, len(sorted_list_EO)):
            hename = sorted_list_EO[i]
            hename = str(hename[:-7])
            if hename in names:
                indices = [i for i, s in enumerate(names) if hename in s]
                x = (float(TeiQueSF_emotionality[int(indices[0])]))
                if x<5.0:
                    x_corr.append(x)
                    y = (float(e[i]))
                    y_corr.append(y)

        #PLOT ALL POINTS
        resultat = stats.pearsonr(x_corr, y_corr)
        minus_5.append(resultat[1])
        
    return minus_5




def max5_return_pval(names, sorted_list_EO, TeiQueSF_emotionality, list_bands):
    all_info = []
    x_corr = []
    y_corr = []
    max_5 = []
    for e in list_bands:
        for i in range(0, len(sorted_list_EO)):
            hename = sorted_list_EO[i]
            hename = str(hename[:-7])
            if hename in names:
                indices = [i for i, s in enumerate(names) if hename in s]
                x = (float(TeiQueSF_emotionality[int(indices[0])]))
                if x>5.0:
                    x_corr.append(x)
                    y = (float(e[i]))
                    y_corr.append(y)

        #PLOT ALL POINTS
        resultat = stats.pearsonr(x_corr, y_corr)
        max_5.append(resultat[1])


    return max_5


def plotbar_info_pval (names, sorted_list_EO, TeiQueSF_emotionality, list_bands):
    labels = ['Betta', 'Gamma', 'ALpha', 'theta', 'delta']
    minus5 = minus_return_pval(names, sorted_list_EO, TeiQueSF_emotionality, list_bands)
    max5 = max5_return_pval(names, sorted_list_EO, TeiQueSF_emotionality, list_bands)

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, minus5, width, label='-5')
    rects2 = ax.bar(x + width/2, max5, width, label='+5')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Two-tailed p-value')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()


    '''def autolabel(rects):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')

    autolabel(rects1)
    autolabel(rects2)'''

    fig.tight_layout()

    plt.show()



def minus_return_pearson(names, sorted_list_EO, TeiQueSF_emotionality, list_bands):
    all_info = []
    x_corr = []
    y_corr = []
    minus_5 = []
    for e in list_bands:
        for i in range(0, len(sorted_list_EO)):
            hename = sorted_list_EO[i]
            hename = str(hename[:-7])
            if hename in names:
                indices = [i for i, s in enumerate(names) if hename in s]
                x = (float(TeiQueSF_emotionality[int(indices[0])]))
                if x<5.0:
                    x_corr.append(x)
                    y = (float(e[i]))
                    y_corr.append(y)

        #PLOT ALL POINTS
        resultat = stats.pearsonr(x_corr, y_corr)
        minus_5.append(resultat[0])
        
    return minus_5




def max5_return_pearson(names, sorted_list_EO, TeiQueSF_emotionality, list_bands):
    all_info = []
    x_corr = []
    y_corr = []
    max_5 = []
    for e in list_bands:
        for i in range(0, len(sorted_list_EO)):
            hename = sorted_list_EO[i]
            hename = str(hename[:-7])
            if hename in names:
                indices = [i for i, s in enumerate(names) if hename in s]
                x = (float(TeiQueSF_emotionality[int(indices[0])]))
                if x>5.0:
                    x_corr.append(x)
                    y = (float(e[i]))
                    y_corr.append(y)

        #PLOT ALL POINTS
        resultat = stats.pearsonr(x_corr, y_corr)
        max_5.append(resultat[0])


    return max_5


def plotbar_info_pearson (names, sorted_list_EO, TeiQueSF_emotionality, list_bands):
    labels = ['Betta', 'Gamma', 'ALpha', 'theta', 'delta']
    minus5 = minus_return_pearson(names, sorted_list_EO, TeiQueSF_emotionality, list_bands)
    max5 = max5_return_pearson(names, sorted_list_EO, TeiQueSF_emotionality, list_bands)

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, minus5, width, label='-5')
    rects2 = ax.bar(x + width/2, max5, width, label='+5')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Pearson corr.')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()


    '''def autolabel(rects):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')

    autolabel(rects1)
    autolabel(rects2)'''

    fig.tight_layout()

    plt.show()