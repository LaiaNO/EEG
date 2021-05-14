
import matplotlib
import matplotlib.pyplot as plt  # plotting
import numpy as np  # linear algebra
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
import scipy
import scipy.interpolate
import statistics
from scipy import signal
import matplotlib.colors as mcolors
import seaborn as sns
import os
from scipy import stats
import statsmodels.api as sm

# Make the correlation of all the data:
# FUNCTION: Print and return the 0: Pearson’s correlation coefficient, 1: Two-tailed p-value
#          Also Plot all the plots with the ...


def calcul(testNames, patientNames, TEST, pacientData, yo, oy):
    # testNames:
    # patientNames:
    # TEST:
    # pacientData: bandPower each one with the data
    # yo & oy:
    savecor = []  # Save the correlation of the xCorr and yCorr
    for e in range(0, len(pacientData)):  # For each frequency band:
        xCorr = []  # Create variables where will go the TEST value
        yCorr = []  # Create variables where will go the pacientData value
        for i in range(0, len(patientNames[e])):  # For each patient:
            hename = patientNames[e][i]  # Select the name of the patient
            # Select only the number witout the extension
            hename = str(hename[:-7])
            if hename in testNames:  # If this is in the testNames:
                # Get the position of the testNames
                indices = [i for i, s in enumerate(testNames) if hename in s]
                # Get the value of the patient in the test selected
                x = (float(TEST[int(indices[0])]))
                xCorr.append(x)  # Add it in to the xCorr
                # Get the value of the patientdata bandPower
                y = (float(pacientData[e][i]))
                yCorr.append(y)  # Add it in to the yCorr
                plt.scatter(x, y, c='g')  # Add it a point in to the plot.
        # Create the line
        z = np.poly1d(np.polyfit(xCorr, yCorr, 1))
        y_len = np.array(len(yCorr))
        xp = np.linspace(yo, oy, y_len)
        y = z(xp)
        plt.plot(xp, y, c='b')

        # PLOT ALL POINTS
        plt.title('Plot show the correlation between test Selected and Band Power')
        plt.show()
        print(stats.pearsonr(xCorr, yCorr))  # Print the stats.pearsonr
        # 0: Pearson’s correlation coefficient, 1: Two-tailed p-value
        savecor.append(stats.pearsonr(xCorr, yCorr))
    return savecor


# BandPower from 4 in TEST


def return_bandpower(testNames, patientNames, TEST, pacientData, minMax, thTIPE):
    val = []  # Save the statistic.median of the y_corr
    for e in range(0, len(pacientData)):  # For each frequency band:
        y_corr = []  # Create variables where will go the pacientData value
        for i in range(0, len(patientNames[e])):  # For each patient:
            hename = patientNames[e][i]  # Select the name of the patient
            # Select only the number witout the extension
            hename = str(hename[:-7])
            if hename in testNames:  # If this is in the testNames:
                # Get the position of the testNames
                indices = [i for i, s in enumerate(testNames) if hename in s]
                # Get the value of the patient in the test selected
                x = (float(TEST[int(indices[0])]))
                if minMax == 'Max':
                    if x > thTIPE:  # If the patient test value is more than 4.
                        # Get the value of the patientdata bandPower
                        y = (float(pacientData[e][i]))
                        y_corr.append(y)  # Add it in to the yCorr
                if minMax == 'Min':
                    if x <= thTIPE:  # If the patient test value is less than 4.
                        # Get the value of the patientdata bandPower
                        y = (float(pacientData[e][i]))
                        y_corr.append(y)  # Add it in to the yCorr
        # PLOT ALL POINTS
        # Median (middle value) of data in that frequency band.
        resultat = statistics.median(y_corr)
        val.append(resultat)  # Add the

    return val

# standar error of the mean


def return_str(testNames, patientNames, TEST, pacientData, minMax, thTIPE):
    stt = []  # Save the scipy.stats.sem of the y_corr
    for e in range(0, len(pacientData)):  # For each frequency band:
        y_corr = []  # Create variables where will go the pacientData value
        for i in range(0, len(patientNames[e])):  # For each patient:
            hename = patientNames[e][i]  # Select the name of the patient
            # Select only the number witout the extension
            hename = str(hename[:-7])
            if hename in testNames:  # If this is in the testNames:
                # Get the position of the testNames
                indices = [i for i, s in enumerate(testNames) if hename in s]
                # Get the value of the patient in the test selected
                x = (float(TEST[int(indices[0])]))
                if minMax == 'Max':
                    if x > thTIPE:  # If the patient test value is more than 4.
                        # Get the value of the patientdata bandPower
                        y = (float(pacientData[e][i]))
                        y_corr.append(y)  # Add it in to the yCorr
                if minMax == 'Min':
                    if x <= thTIPE:  # If the patient test value is less than 4.
                        # Get the value of the patientdata bandPower
                        y = (float(pacientData[e][i]))
                        y_corr.append(y)  # Add it in to the yCorr

        # The standard error of the mean in the sample(s).
        resultat = scipy.stats.sem(y_corr)
        stt.append(resultat)  # Add the standard error of the mean
    return stt


def plotbar_info_bandpower1(names, sorted_list_EO, Test, list_bands, thTIPE):
    # chanels1_delta_EO, chanels1_theta_EO, chanels1_alpha_EO, chanels1_beta_EO, chanels1_gama_l_EO, chanels1_gama_u_EO
    labels = ['Delta', 'Theta', 'Alpha', 'Beta', 'Gama Lower']
    # the median of the minus 4
    minus5 = return_bandpower(names, sorted_list_EO,
                              Test, list_bands, 'Min', thTIPE)
    max5 = return_bandpower(names, sorted_list_EO, Test,
                            list_bands, 'Max', thTIPE)  # the median of the max 4

    strmax = return_str(names, sorted_list_EO, Test, list_bands,
                        'Max', thTIPE)  # standar error of the mean
    strmin = return_str(names, sorted_list_EO, Test, list_bands,
                        'Min', thTIPE)  # standar error of the mean

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    ax.bar(x - width/2, minus5, width, yerr=strmin,
           capsize=0.05, label='-'+str(thTIPE))
    ax.bar(x + width/2, max5, width, yerr=strmax,
           capsize=0.05, label='+'+str(thTIPE))

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Band Power')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    fig.tight_layout()
    plt.show()

    return max5


def return_rac(testNames, patientNames, TEST, pacientData, minMax, thTIPE):
    info_chan = []  # Save the y_corr < 4 in TEST
    for e in range(0, len(pacientData)):  # For each frequency band:
        y_corr = []  # Create variables where will go the pacientData value
        for i in range(0, len(patientNames[e])):  # For each patient:
            hename = patientNames[e][i]  # Select the name of the patient
            # Select only the number witout the extension
            hename = str(hename[:-7])
            if hename in testNames:  # If this is in the testNames:
                # Get the position of the testNames
                indices = [i for i, s in enumerate(testNames) if hename in s]
                # Get the value of the patient in the test selected
                x = (float(TEST[int(indices[0])]))
                if minMax == 'Max':
                    if x > thTIPE:  # If the patient test value is more than 4.
                        # Get the value of the patientdata bandPower
                        y = (float(pacientData[e][i]))
                        y_corr.append(y)  # Add it in to the yCorr
                if minMax == 'Min':
                    if x < thTIPE:  # If the patient test value is less than 4.
                        # Get the value of the patientdata bandPower
                        y = (float(pacientData[e][i]))
                        y_corr.append(y)  # Add it in to the yCorr
        info_chan.append(y_corr)  # Save it
    return info_chan


def plotbar_info_rac(names, sorted_list_EO, TeiQueSF_emotionality, list_bands, thTIPE):
    # chanels1_delta_EO, chanels1_theta_EO, chanels1_alpha_EO, chanels1_beta_EO, chanels1_gama_l_EO, chanels1_gama_u_EO
    # labels = ['Delta', 'Theta', 'Alpha', 'Beta', 'Gama Lower']

    minus4 = return_rac(names, sorted_list_EO,
                        TeiQueSF_emotionality, list_bands, 'Min', thTIPE)
    max4 = return_rac(names, sorted_list_EO,
                      TeiQueSF_emotionality, list_bands, 'Max', thTIPE)

    resu = []
    for i in range(0, len(max4)):
        # The test statistic under the large-sample approximation that the rank sum statistic is normally distributed; pvaluefloat The two-sided p-value of the test.
        p = scipy.stats.ranksums(minus4[i], max4[i])
        print(p)
        resu.append(p)
    return resu
