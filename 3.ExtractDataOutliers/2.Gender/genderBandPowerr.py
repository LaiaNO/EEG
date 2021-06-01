
import matplotlib
import matplotlib.pyplot as plt  # plotting
import numpy as np  # linear algebra
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
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

# Make the correlation of all the data:
# FUNCTION: Print and return the 0: Pearson’s correlation coefficient, 1: Two-tailed p-value
#          Of all the women and all the men


def calcul(testNames, patientNames, TEST, gender, pacientData, yo, oy):
    # testNames:
    # patientNames:
    # TEST:
    # Gender:
    # pacientData: bandPower each one with the data
    # yo & oy:
    savecor_f = []  # Save the correlation of the xCorr and yCorr for Female
    savecor_m = []  # Save the correlation of the xCorr and yCorr for Male
    for e in range(0, len(pacientData)):  # For each frequency band:
        xCorr_m = []  # Create variables where will go the TEST value for Male
        yCorr_m = []  # Create variables where will go the pacientData value for Male
        xCorr_f = []  # Create variables where will go the TEST value for Female
        yCorr_f = []  # Create variables where will go the pacientData value for Female
        for i in range(0, len(patientNames[e])):  # For each patient:
            hename = patientNames[e][i]  # Select the name of the patient
            # Select only the number witout the extension
            hename = str(hename[:-7])
            if hename in testNames:  # If this is in the testNames:
                # Get the position of the testNames
                indices = [i for i, s in enumerate(testNames) if hename in s]
                # Get the value of the patient in the test selected
                x = (float(TEST[int(indices[0])]))
                # Get the gender of the patient
                g = (float(gender[int(indices[0])]))
                if str(g) == '1.0':
                    xCorr_f.append(x)
                    y = (float(pacientData[e][i]))
                    yCorr_f.append(y)
                if str(g) == '2.0':
                    xCorr_m.append(x)
                    y = (float(pacientData[e][i]))
                    yCorr_m.append(y)

        print(stats.pearsonr(xCorr_f, yCorr_f))  # Print the stats.pearsonr
        # 0: Pearson’s correlation coefficient, 1: Two-tailed p-value
        savecor_f.append(stats.pearsonr(xCorr_f, yCorr_f))

        print(stats.pearsonr(xCorr_m, yCorr_m))  # Print the stats.pearsonr
        # 0: Pearson’s correlation coefficient, 1: Two-tailed p-value
        savecor_m.append(stats.pearsonr(xCorr_m, yCorr_m))
    return savecor_f, savecor_m

# BandPower Minus from 4 in TEST


def return_bandpower1(testNames, patientNames, TEST, gender, pacientData, minMax):
    minusF = []  # Save the statistic.median of the y_corr
    minusM = []  # Save the statistic.median of the y_corr
    for e in range(0, len(pacientData)):  # For each frequency band:
        y_corr_f = []  # Create variables where will go the pacientData value
        y_corr_m = []  # Create variables where will go the pacientData value
        for i in range(0, len(patientNames[e])):  # For each patient:
            hename = patientNames[e][i]  # Select the name of the patient
            # Select only the number witout the extension
            hename = str(hename[:-7])
            if hename in testNames:  # If this is in the testNames:
                # Get the position of the testNames
                indices = [i for i, s in enumerate(testNames) if hename in s]
                # Get the value of the patient in the test selected
                x = (float(TEST[int(indices[0])]))
                # Get the gender of the patient
                g = (float(gender[int(indices[0])]))
                if minMax == 'Max':
                    if str(g) == '1.0' and x > 4.5:
                        y = (float(pacientData[e][i]))
                        y_corr_f.append(y)
                    if str(g) == '2.0' and x > 4.5:
                        y = (float(pacientData[e][i]))
                        y_corr_m.append(y)
                if minMax == 'Min':
                    if str(g) == '1.0' and x < 4.5:
                        y = (float(pacientData[e][i]))
                        y_corr_f.append(y)
                    if str(g) == '2.0' and x < 4.5:
                        y = (float(pacientData[e][i]))
                        y_corr_m.append(y)

        # PLOT ALL POINTS
        resultatF = statistics.median(y_corr_f)
        minusF.append(resultatF)

        resultatM = statistics.median(y_corr_m)
        minusM.append(resultatM)

    return minusF, minusM

# standar error of the mean -> np.std(data, ddof=1) / np.sqrt(np.size(data))


def return_str1(testNames, patientNames, TEST, gender, pacientData, minMax):
    stt_f = []
    stt_m = []
    for e in range(0, len(pacientData)):  # For each frequency band:
        y_corr_f = []  # Create variables where will go the pacientData value
        y_corr_m = []  # Create variables where will go the pacientData value
        for i in range(0, len(patientNames[e])):  # For each patient:
            hename = patientNames[e][i]  # Select the name of the patient
            # Select only the number witout the extension
            hename = str(hename[:-7])
            if hename in testNames:  # If this is in the testNames:
                # Get the position of the testNames
                indices = [i for i, s in enumerate(testNames) if hename in s]
                # Get the value of the patient in the test selected
                x = (float(TEST[int(indices[0])]))
                # Get the gender of the patient
                g = (float(gender[int(indices[0])]))
                if minMax == 'Max':
                    if str(g) == '1.0' and x > 4.5:
                        y = (float(pacientData[e][i]))
                        y_corr_f.append(y)
                    if str(g) == '2.0' and x > 4.5:
                        y = (float(pacientData[e][i]))
                        y_corr_m.append(y)
                if minMax == 'Min':
                    if str(g) == '1.0' and x < 4.5:
                        y = (float(pacientData[e][i]))
                        y_corr_f.append(y)
                    if str(g) == '2.0' and x < 4.5:
                        y = (float(pacientData[e][i]))
                        y_corr_m.append(y)

        # PLOT ALL POINTS
        resultat_m = scipy.stats.sem(y_corr_m)
        stt_m.append(resultat_m)

        resultat_f = scipy.stats.sem(y_corr_f)
        stt_f.append(resultat_f)

    return stt_f, stt_m


def plotbar_info_bandpower1(names, sorted_list_EO, TeiQueSF_emotionality, gender, list_bands):
    #chanels1_delta_EO, chanels1_theta_EO, chanels1_alpha_EO, chanels1_beta_EO, chanels1_gama_l_EO, chanels1_gama_u_EO
    labels = ['Delta', 'Theta', 'Alpha', 'Beta', 'Gama Lower']
    minus5_f, minus5_m = return_bandpower1(
        names, sorted_list_EO, TeiQueSF_emotionality, gender, list_bands, 'Min')
    max_5_f, max_5_m = return_bandpower1(
        names, sorted_list_EO, TeiQueSF_emotionality, gender, list_bands, 'Max')
    # print(minus5)
    strfMax, strmMax = return_str1(
        names, sorted_list_EO, TeiQueSF_emotionality, gender, list_bands, 'Max')
    strfMin, strmMin = return_str1(
        names, sorted_list_EO, TeiQueSF_emotionality, gender, list_bands, 'Min')

    x = np.arange(len(labels))  # the label locations
    width = 0.20  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x-0.2, minus5_f, width, yerr=strfMin,
                    capsize=0.05, label='-4 women')
    rects3 = ax.bar(x, max_5_f, width, yerr=strfMax,
                    capsize=0.05, label='+4 women')

    rects2 = ax.bar(x+0.2, minus5_m, width, yerr=strmMin,
                    capsize=0.05, label='-4 men')
    rects4 = ax.bar(x+0.4, max_5_m, width, yerr=strmMax,
                    capsize=0.05, label='+4 men')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Band Power')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    fig.tight_layout()

    plt.show()

    return max_5_f, max_5_m


def return_rac(testNames, patientNames, TEST, gender, pacientData, minMax):
    info_chan_f = []
    info_chan_m = []
    for e in range(0, len(pacientData)):  # For each frequency band:
        y_corr_f = []  # Create variables where will go the pacientData value
        y_corr_m = []  # Create variables where will go the pacientData value
        for i in range(0, len(patientNames[e])):  # For each patient:
            hename = patientNames[e][i]  # Select the name of the patient
            # Select only the number witout the extension
            hename = str(hename[:-7])
            if hename in testNames:  # If this is in the testNames:
                # Get the position of the testNames
                indices = [i for i, s in enumerate(testNames) if hename in s]
                # Get the value of the patient in the test selected
                x = (float(TEST[int(indices[0])]))
                # Get the gender of the patient
                g = (float(gender[int(indices[0])]))
                if minMax == 'Max':
                    if str(g) == '1.0' and x > 4.5:
                        y = (float(pacientData[e][i]))
                        y_corr_f.append(y)
                    if str(g) == '2.0' and x > 4.5:
                        y = (float(pacientData[e][i]))
                        y_corr_m.append(y)
                if minMax == 'Min':
                    if str(g) == '1.0' and x < 4.5:
                        y = (float(pacientData[e][i]))
                        y_corr_f.append(y)
                    if str(g) == '2.0' and x < 4.5:
                        y = (float(pacientData[e][i]))
                        y_corr_m.append(y)
        info_chan_m.append(y_corr_m)
        info_chan_f.append(y_corr_f)

    return info_chan_f, info_chan_m


def plotbar_info_rac(names, sorted_list_EO, TeiQueSF_emotionality, gender, list_bands):
    #chanels1_delta_EO, chanels1_theta_EO, chanels1_alpha_EO, chanels1_beta_EO, chanels1_gama_l_EO, chanels1_gama_u_EO
    labels = ['Delta', 'Theta', 'Alpha', 'Beta', 'Gama Lower']

    minusF, minusM = return_rac(
        names, sorted_list_EO, TeiQueSF_emotionality, gender, list_bands, 'Min')
    maxF, maxM = return_rac(names, sorted_list_EO,
                            TeiQueSF_emotionality, gender, list_bands, 'Max')

    resuM = []
    for i in range(0, len(maxM)):
        p = scipy.stats.ranksums(minusM[i], maxM[i])
        resuM.append(p)
    print('Men')
    print(resuM)

    resuF = []
    for i in range(0, len(maxF)):
        p = scipy.stats.ranksums(minusF[i], maxF[i])
        resuF.append(p)
    print('Women')
    print(resuF)

    return resuF, resuM
