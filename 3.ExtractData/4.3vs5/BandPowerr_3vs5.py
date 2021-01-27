
import matplotlib
import matplotlib.pyplot as plt # plotting
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
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

#Make the correlation of all the data:
def calcul(testNames, patientNames, TEST, pacientData, yo, oy):
    #testNames: 
    #patientNames:
    #TEST:
    #pacientData: bandPower each one with the data
    #yo & oy: 
    savecor = [] #Save the correlation of the xCorr and yCorr
    for e in range(0,len(pacientData)):       #For each frequency band:
        xCorr = []   #Create variables where will go the TEST value
        yCorr = []   #Create variables where will go the pacientData value
        for i in range(0, len(patientNames)): #For each patient:
            hename = patientNames[i]   #Select the name of the patient
            hename = str(hename[:-7])  #Select only the number witout the extension 
            if hename in testNames:    #If this is in the testNames:
                indices = [i for i, s in enumerate(testNames) if hename in s] #Get the position of the testNames
                x = (float(TEST[int(indices[0])])) #Get the value of the patient in the test selected
                xCorr.append(x) #Add it in to the xCorr
                y = (float(pacientData[e][i]))     #Get the value of the patientdata bandPower
                yCorr.append(y) #Add it in to the yCorr
                plt.scatter(x,y, c='g') #Add it a point in to the plot.
        #Create the line
        z = np.poly1d(np.polyfit(xCorr, yCorr, 1))
        y_len = np.array(len(yCorr))
        xp = np.linspace(yo, oy, y_len)
        y = z(xp)
        plt.plot(xp, y, c='b')

        #PLOT ALL POINTS
        plt.title('Plot show the correlation between test Selected and Band Power')
        plt.show()
        print(stats.pearsonr(xCorr, yCorr)) #Print the stats.pearsonr
        savecor.append(stats.pearsonr(xCorr, yCorr)) #0: Pearson’s correlation coefficient, 1: Two-tailed p-value
    return savecor

def calcul_More4(testNames, patientNames, TEST, pacientData):
    #testNames: 
    #patientNames:
    #TEST:
    #pacientData: bandPower each one with the data
    #yo & oy: 
    savecor = [] #Save the correlation of the xCorr and yCorr
    for e in range(0,len(pacientData)):       #For each frequency band:
        xCorr = []   #Create variables where will go the TEST value
        yCorr = []   #Create variables where will go the pacientData value
        for i in range(0, len(patientNames)): #For each patient:
            hename = patientNames[i]   #Select the name of the patient
            hename = str(hename[:-7])  #Select only the number witout the extension 
            if hename in testNames:    #If this is in the testNames:
                indices = [i for i, s in enumerate(testNames) if hename in s] #Get the position of the testNames
                x = (float(TEST[int(indices[0])])) #Get the value of the patient in the test selected
                if x>4.0:
                    xCorr.append(x) #Add it in to the xCorr
                    y = (float(pacientData[e][i]))     #Get the value of the patientdata bandPower
                    yCorr.append(y) #Add it in to the yCorr
        savecor.append(stats.pearsonr(xCorr, yCorr)) #0: Pearson’s correlation coefficient, 1: Two-tailed p-value
    return savecor

#BandPower Minus from 4 in TEST
def minus_return_bandpower(testNames, patientNames, TEST, pacientData):
    minusVal = [] #Save the statistic.median of the y_corr
    for e in pacientData: #For each frequency band:
        y_corr = []       #Create variables where will go the pacientData value
        for i in range(0, len(patientNames)): #For each patient:
            hename = patientNames[i]  #Select the name of the patient
            hename = str(hename[:-7]) #Select only the number witout the extension
            if hename in testNames:   #If this is in the testNames:
                indices = [i for i, s in enumerate(testNames) if hename in s] #Get the position of the testNames
                x = (float(TEST[int(indices[0])])) #Get the value of the patient in the test selected
                if x<4.0: #If the patient test value is less than 4.
                    y = (float(e[i])) #Get the value of the patientdata bandPower
                    y_corr.append(y)  #Add it in to the yCorr

        #PLOT ALL POINTS
        resultat = statistics.median(y_corr) #Median (middle value) of data in that frequency band.
        minusVal.append(resultat) #Add the median
        
    return minusVal

#BandPower More from 4 in TEST
def max_return_bandpower(testNames, patientNames, TEST, pacientData):
    maxVal = [] #Save the statistic.median of the y_corr
    for e in pacientData: #For each frequency band:
        y_corr = []       #Create variables where will go the pacientData value
        for i in range(0, len(patientNames)): #For each patient:
            hename = patientNames[i]  #Select the name of the patient
            hename = str(hename[:-7]) #Select only the number witout the extension
            if hename in testNames:   #If this is in the testNames:
                indices = [i for i, s in enumerate(testNames) if hename in s] #Get the position of the testNames
                x = (float(TEST[int(indices[0])])) #Get the value of the patient in the test selected
                if x>4.0: #If the patient test value is more than 4.
                    y = (float(e[i])) #Get the value of the patientdata bandPower
                    y_corr.append(y)  #Add it in to the yCorr

        #PLOT ALL POINTS
        resultat = statistics.median(y_corr) #Median (middle value) of data in that frequency band.
        maxVal.append(resultat) #Add the median
        
    return maxVal

#standar error of the mean 
def max_return_str(testNames, patientNames, TEST, pacientData):
    stt = [] #Save the scipy.stats.sem of the y_corr
    for e in pacientData: #For each frequency band:
        y_corr = []       #Create variables where will go the pacientData value
        for i in range(0, len(patientNames)): #For each patient:
            hename = patientNames[i]  #Select the name of the patient
            hename = str(hename[:-7]) #Select only the number witout the extension
            if hename in testNames:   #If this is in the testNames:
                indices = [i for i, s in enumerate(testNames) if hename in s] #Get the position of the testNames
                x = (float(TEST[int(indices[0])])) #Get the value of the patient in the test selected
                if x>4.0: #If the patient test value is more than 4.
                    y = (float(e[i])) #Get the value of the patientdata bandPower
                    y_corr.append(y)  #Add it in to the yCorr

        resultat = scipy.stats.sem(y_corr) #The standard error of the mean in the sample(s).
        stt.append(resultat) #Add the standard error of the mean
    return stt

#standar error of the mean 
def min_return_str1(testNames, patientNames, TEST, pacientData):
    stt = [] #Save the scipy.stats.sem of the y_corr
    for e in pacientData: #For each frequency band:
        y_corr = []       #Create variables where will go the pacientData value
        for i in range(0, len(patientNames)): #For each patient:
            hename = patientNames[i]  #Select the name of the patient
            hename = str(hename[:-7]) #Select only the number witout the extension
            if hename in testNames:   #If this is in the testNames:
                indices = [i for i, s in enumerate(testNames) if hename in s] #Get the position of the testNames
                x = (float(TEST[int(indices[0])])) #Get the value of the patient in the test selected
                if x<4.0: #If the patient test value is more than 4.
                    y = (float(e[i])) #Get the value of the patientdata bandPower
                    y_corr.append(y)  #Add it in to the yCorr

        resultat = scipy.stats.sem(y_corr) #f = np.std(y_corr, ddof=1) / np.sqrt(np.size(y_corr))
        '''y_mean = np.mean(y_corr)
        y_std = np.std(y_corr)
        meen.append(y_mean)'''
        stt.append(resultat) #Save
    return stt

def plotbar_info_bandpower1 (names, sorted_list_EO, TeiQueSF_emotionality, list_bands):
    #chanels1_delta_EO, chanels1_theta_EO, chanels1_alpha_EO, chanels1_beta_EO, chanels1_gama_l_EO, chanels1_gama_u_EO
    labels = ['Delta', 'Theta', 'Alpha', 'Beta', 'Gama Lower']
    minus5 = minus_return_bandpower(names, sorted_list_EO, TeiQueSF_emotionality, list_bands) #the median of the minus 4
    max5 = max_return_bandpower(names, sorted_list_EO, TeiQueSF_emotionality, list_bands) #the median of the max 4
    str1 = max_return_str(names, sorted_list_EO, TeiQueSF_emotionality, list_bands) #standar error of the mean 
    str2 = min_return_str1(names, sorted_list_EO, TeiQueSF_emotionality, list_bands)#standar error of the mean 

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars


    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, minus5, width, yerr=str2,capsize=0.05,label='-3')
    rects2 = ax.bar(x + width/2, max5, width, yerr=str1,capsize=0.05, label='+5')


    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Band Power')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    fig.tight_layout()
    plt.show()

    return max5


def minus_return_rac(testNames, patientNames, TEST, pacientData):
    info_chan = [] #Save the y_corr < 4 in TEST
    for e in pacientData: #For each frequency band:
        y_corr = [] #Create variables where will go the pacientData value
        for i in range(0, len(patientNames)): #For each patient:
            hename = patientNames[i]  #Select the name of the patient
            hename = str(hename[:-7]) #Select only the number witout the extension
            if hename in testNames:   #If this is in the testNames:
                indices = [i for i, s in enumerate(testNames) if hename in s] #Get the position of the testNames
                x = (float(TEST[int(indices[0])])) #Get the value of the patient in the test selected
                if x<4.0: #If the patient test value is more than 4.
                    y = (float(e[i])) #Get the value of the patientdata bandPower
                    y_corr.append(y)  #Add it in to the yCorr
        info_chan.append(y_corr)  #Save it   
    return info_chan




def max_return_rac(testNames, patientNames, TEST, pacientData):
    info_chan = [] #Save the y_corr < 4 in TEST
    for e in pacientData: #For each frequency band:
        y_corr = [] #Create variables where will go the pacientData value
        for i in range(0, len(patientNames)): #For each patient:
            hename = patientNames[i]  #Select the name of the patient
            hename = str(hename[:-7]) #Select only the number witout the extension
            if hename in testNames:   #If this is in the testNames:
                indices = [i for i, s in enumerate(testNames) if hename in s] #Get the position of the testNames
                x = (float(TEST[int(indices[0])])) #Get the value of the patient in the test selected
                if x>4.0: #If the patient test value is more than 4.
                    y = (float(e[i])) #Get the value of the patientdata bandPower
                    y_corr.append(y)  #Add it in to the yCorr
        info_chan.append(y_corr)  #Save it   
    return info_chan


def plotbar_info_rac (names, sorted_list_EO, TeiQueSF_emotionality, list_bands):
    #chanels1_delta_EO, chanels1_theta_EO, chanels1_alpha_EO, chanels1_beta_EO, chanels1_gama_l_EO, chanels1_gama_u_EO
    labels = ['Delta', 'Theta', 'Alpha', 'Beta', 'Gama Lower']

    minus4 = minus_return_rac(names, sorted_list_EO, TeiQueSF_emotionality, list_bands)
    max4 = max_return_rac(names, sorted_list_EO, TeiQueSF_emotionality, list_bands)

    resu = []
    for i in range(0,len(max4)):
        p = scipy.stats.ranksums(minus4[i],max4[i]) #The test statistic under the large-sample approximation that the rank sum statistic is normally distributed; pvaluefloat The two-sided p-value of the test.
        print (p)
        resu.append(p)
    return resu