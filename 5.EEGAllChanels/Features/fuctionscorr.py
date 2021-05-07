
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

# Make the correlation of all the data:
# FUNCTION: Print and return the 0: Pearson’s correlation coefficient, 1: Two-tailed p-value
#          Also Plot all the plots with the ...


def calcul(testNames, Names, TEST, DataEEG, yo, oy):
    # testNames:
    # patientNames:
    # TEST:
    # pacientData: bandPower each one with the data
    # yo & oy:

    savecor = []  # Save the correlation of the xCorr and yCorr
    for r in range(0, 3):
        xCorr = []
        yCorr = []
        for i in range(0, len(Names)):  # For each patient:
            hename = Names[i]  # Select the name of the patient
            if hename in testNames:  # If this is in the testNames:
                # Get the position of the testNames
                indices = [i for i, s in enumerate(testNames) if hename in s]
                # Get the value of the patient in the test selected
                x = (float(TEST[int(indices[0])]))
                xCorr.append(x)  # Add it in to the xCorr
                # Get the value of the patientdata bandPower
                y = (float(DataEEG[i][r]))
                yCorr.append(y)  # Add it in to the yCorr
                plt.scatter(x, y, c='g')  # Add it a point in to the plot.
        if r == 0:
            print('Valence')
        if r == 1:
            print('Arousal')
        if r == 2:
            print('Motivation')

        z = np.poly1d(np.polyfit(xCorr, yCorr, 1))
        y_len = np.array(len(yCorr))
        xp = np.linspace(yo, oy, y_len)
        y = z(xp)
        plt.plot(xp, y, c='b')
        # PLOT ALL POINTS
        plt.title('Plot show the correlation between test Selected and -')
        plt.show()

        print(stats.pearsonr(xCorr, yCorr))  # Print the stats.pearsonr
        # 0: Pearson’s correlation coefficient, 1: Two-tailed p-value
        savecor.append(stats.pearsonr(xCorr, yCorr))
    return savecor


# BandPower from 4 in TEST
def return_features(testNames, Names, TEST, DataEEG, minMax, val, thTIPE):
    y_corr = []
    for i in range(0, len(Names)):  # For each patient:
        hename = Names[i]
        if hename in testNames:  # If this is in the testNames:
            # Get the position of the testNames
            indices = [i for i, s in enumerate(testNames) if hename in s]
            # Get the value of the patient in the test selected
            x = (float(TEST[int(indices[0])]))
            if minMax == 'Max':
                if x >= thTIPE:  # If the patient test value is more than 4.
                    # Get the value of the patientdata bandPower
                    y = (float(DataEEG[i][val]))
                    y_corr.append(y)  # Add it in to the yCorr
            if minMax == 'Min':
                if x < thTIPE:  # If the patient test value is less than 4.
                    # Get the value of the patientdata bandPower
                    y = (float(DataEEG[i][val]))
                    y_corr.append(y)  # Add it in to the yCorr

    # PLOT ALL POINTS
    # Median (middle value) of data in that frequency band.
    resultat = statistics.median(y_corr)

    return resultat

# standar error of the mean


def return_str(testNames, Names, TEST, DataEEG, minMax, val, th):
    stt = []  # Save the scipy.stats.sem of the y_corr
    y_corr = []  # Create variables where will go the pacientData value
    for i in range(0, len(Names)):  # For each patient:
        hename = Names[i]  # Select the name of the patient
        if hename in testNames:  # If this is in the testNames:
            # Get the position of the testNames
            indices = [i for i, s in enumerate(testNames) if hename in s]
            # Get the value of the patient in the test selected
            x = (float(TEST[int(indices[0])]))
            if minMax == 'Max':
                if x >= th:  # If the patient test value is more than 4.
                    # Get the value of the patientdata bandPower
                    y = (float(DataEEG[i][val]))
                    y_corr.append(y)  # Add it in to the yCorr
            if minMax == 'Min':
                if x < th:  # If the patient test value is less than 4.
                    # Get the value of the patientdata bandPower
                    y = (float(DataEEG[i][val]))
                    y_corr.append(y)  # Add it in to the yCorr

    # The standard error of the mean in the sample(s).
    resultat = scipy.stats.sem(y_corr)
    stt.append(resultat)  # Add the standard error of the mean
    return stt


def plotbar_info(testNames, Names, TEST, DataEEG, thTIPE):
    data = []
    #chanels1_delta_EO, chanels1_theta_EO, chanels1_alpha_EO, chanels1_beta_EO, chanels1_gama_l_EO, chanels1_gama_u_EO
    labels = ['Valence']
    print('Valence')
    # the median of the minus 4
    minus5 = return_features(testNames, Names, TEST, DataEEG, 'Min', 0, thTIPE)
    max5 = return_features(testNames, Names, TEST, DataEEG,
                           'Max', 0, thTIPE)  # the median of the max 4

    str1 = return_str(testNames, Names, TEST, DataEEG,
                      'Min', 0, thTIPE)  # standar error of the mean
    str2 = return_str(testNames, Names, TEST, DataEEG,
                      'Max', 0, thTIPE)  # standar error of the mean

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, minus5, width,
                    yerr=str1, capsize=0.05, label='-4')
    rects2 = ax.bar(x + width/2, max5, width,
                    yerr=str2, capsize=0.05, label='+4')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Band Power')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    fig.tight_layout()
    plt.show()
    data.append(max5)

    # 1 feature
    labels = ['Arousal']
    print('Arousal')
    # the median of the minus 4
    minus5 = return_features(testNames, Names, TEST, DataEEG, 'Min', 1, thTIPE)
    max5 = return_features(testNames, Names, TEST, DataEEG,
                           'Max', 1, thTIPE)  # the median of the max 4

    str1 = return_str(testNames, Names, TEST, DataEEG,
                      'Min', 1, thTIPE)  # standar error of the mean
    str2 = return_str(testNames, Names, TEST, DataEEG,
                      'Max', 1, thTIPE)  # standar error of the mean

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, minus5, width,
                    yerr=str1, capsize=0.05, label='-4')
    rects2 = ax.bar(x + width/2, max5, width,
                    yerr=str2, capsize=0.05, label='+4')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Band Power')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    fig.tight_layout()
    plt.show()
    data.append(max5)

    # feature 2
    labels = ['Motivation']
    print('Motivation')
    # the median of the minus 4
    minus5 = return_features(testNames, Names, TEST, DataEEG, 'Min', 2, thTIPE)
    max5 = return_features(testNames, Names, TEST, DataEEG,
                           'Max', 2, thTIPE)  # the median of the max 4

    str1 = return_str(testNames, Names, TEST, DataEEG,
                      'Min', 2, thTIPE)  # standar error of the mean
    str2 = return_str(testNames, Names, TEST, DataEEG,
                      'Max', 2, thTIPE)  # standar error of the mean

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, minus5, width,
                    yerr=str1, capsize=0.05, label='-4')
    rects2 = ax.bar(x + width/2, max5, width,
                    yerr=str2, capsize=0.05, label='+4')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Band Power')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    fig.tight_layout()
    plt.show()
    data.append(max5)

    return data


def return_rac(testNames, Names, TEST, DataEEG, minMax, val, thTIPE):
    info_chan = []  # Save the y_corr < 4 in TEST
    y_corr = []  # Create variables where will go the pacientData value

    for i in range(0, len(Names)):  # For each patient:
        hename = Names[i]  # Select the name of the patient
        if hename in testNames:  # If this is in the testNames:
            # Get the position of the testNames
            indices = [i for i, s in enumerate(testNames) if hename in s]
            # Get the value of the patient in the test selected
            x = (float(TEST[int(indices[0])]))
            if minMax == 'Max':
                if x >= thTIPE:  # If the patient test value is more than 4.
                    # Get the value of the patientdata bandPower
                    y = (float(DataEEG[i][val]))
                    y_corr.append(y)  # Add it in to the yCorr
            if minMax == 'Min':
                if x < thTIPE:  # If the patient test value is less than 4.
                    # Get the value of the patientdata bandPower
                    y = (float(DataEEG[i][val]))
                    y_corr.append(y)  # Add it in to the yCorr
    info_chan.append(y_corr)  # Save it
    return info_chan


def plotbar_info_rac(testNames, Names, TEST, DataEEG, thTIPE):
    #chanels1_delta_EO, chanels1_theta_EO, chanels1_alpha_EO, chanels1_beta_EO, chanels1_gama_l_EO, chanels1_gama_u_EO
    resu = []
    print('Valence')
    minus4 = return_rac(testNames, Names, TEST, DataEEG, 'Min', 0, thTIPE)
    max4 = return_rac(testNames, Names, TEST, DataEEG, 'Max', 0, thTIPE)
    for i in range(0, len(max4)):
        # The test statistic under the large-sample approximation that the rank sum statistic is normally distributed; pvaluefloat The two-sided p-value of the test.
        p = scipy.stats.ranksums(minus4[i], max4[i])
        print(p)
        resu.append(p)

    print('Arousal')
    minus4 = return_rac(testNames, Names, TEST, DataEEG, 'Min', 1, thTIPE)
    max4 = return_rac(testNames, Names, TEST, DataEEG, 'Max', 1, thTIPE)
    for i in range(0, len(max4)):
        # The test statistic under the large-sample approximation that the rank sum statistic is normally distributed; pvaluefloat The two-sided p-value of the test.
        p = scipy.stats.ranksums(minus4[i], max4[i])
        print(p)
        resu.append(p)

    print('Motivation')
    minus4 = return_rac(testNames, Names, TEST, DataEEG, 'Min', 2, thTIPE)
    max4 = return_rac(testNames, Names, TEST, DataEEG, 'Max', 2, thTIPE)
    for i in range(0, len(max4)):
        # The test statistic under the large-sample approximation that the rank sum statistic is normally distributed; pvaluefloat The two-sided p-value of the test.
        p = scipy.stats.ranksums(minus4[i], max4[i])
        print(p)
        resu.append(p)

    return resu
