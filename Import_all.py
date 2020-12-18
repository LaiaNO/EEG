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

from Definicions import upperchanel
from Definicions import group_inf
from Definicions import mediumchanels
from Definicions import opteciogrups
from Definicions import find_nearest
from Definicions import grabt

def Import_Patients(basePATH):

    ### DECLARE THE 12 CHANELS GROUPS WE WANT TO REDUCE THE 62 CHANELS
    group1= ["FP2", "AFZ", "FZ", "FCZ", "CZ"]
    group2= ["F7", "F5", "F3", "F1", "AF7", "AF3", "FP1"]
    group3= ["FP2", "AF4", "AF8", "F2", "F4", "F6", "F8"]
    group4= ["FT7", "T7", "TP7"]
    group5= ["FC5", "FC3", "FC1", "C5", "C3", "C1"]
    group6= ["CP5", "CP3", "CP1", "P7", "P5", "P3", "P1"]
    group7= ["PO7", "PO3", "O1", "PO9"]
    group8= ["PO4", "PO8", "O2", "PO10"]
    group9= ["P2", "P4", "P6", "P8", "CP2", "CP4", "CP6"]
    group10= ["FT8", "T8", "TP8"]
    group11= ["CPZ", "PZ", "POZ", "OZ", "IZ"]
    group12= ["FC2", "FC4", "FC6", "C2", "C4", "C6"]




    ###SAVE THIS 12 IN A MATRIX TO BE ITERATED LATER
    groups = [group1, group2, group3, group4, group5, group6, group7, group8, group9, group10, group11, group12]

    ###DECLARE WHAT TIPE ARE THIS CHANELS TO GO FROM NUMPY TO MNE
    titlegrup = ["ch1","ch2","ch3","ch4","ch5","ch6","ch7","ch8","ch9","ch10","ch11","ch12"]
    ct_ty = ["eeg", "eeg", "eeg", "eeg", "eeg", "eeg", "eeg", "eeg", "eeg", "eeg", "eeg", "eeg"]
    sfreq = 250
    info = mne.create_info(ch_names=titlegrup, sfreq=sfreq, ch_types=ct_ty)





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

    for i in sorted_list_EO:
        #Load the doc
        x=mne.io.read_raw_eeglab(basePATH+i, preload=True, verbose=True)
        #GET DATA
        data = x._data
        datafinal = data[:119200]
        #GET CHANELS
        chanles_names = x.ch_names
        #REDUCE CHANELS TO 12
        EO = opteciogrups(chanles_names, groups, datafinal)

        #Select the same but with EC
        EC_name = [sorted_list_EC.index(e) for e in sorted_list_EC if i[:-6] in e]
        EC_name_PATH = sorted_list_EC[EC_name[0]]
        #LOAD EC
        x2=mne.io.read_raw_eeglab(basePATH+EC_name_PATH, preload=True, verbose=True)
        #GET DATA
        data2 = x2._data
        datafinal2 = data2[:119200]
        #GET CHANELS
        chanles_names2 = x2.ch_names
        #REDUCE CHANELS TO 12
        EC = opteciogrups(chanles_names2, groups, datafinal2)

        #save
        EO_EC_P = []
        EO_EC_P.append(EO)
        EO_EC_P.append(EC)
        EO_EC_Pacients.append(EO_EC_P)
    
    return(sorted_list_EO, sorted_list_EC, EO_EC_Pacients)
