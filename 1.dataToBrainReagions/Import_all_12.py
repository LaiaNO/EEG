import numpy as np  # linear algebra
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
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


'''FOR THE INITIAL DATA TREATMENT'''


# DEFINITION UPPER CHANNEL (to have all at the same type):
def upperchanel(chanlesnames):
    yourList = [x.upper() for x in chanlesnames]  # Upper the name chanels
    return yourList


# DEFINITION TO FIND IF THERE ARE ALL THE CHANELS WE DECLARE IN THE 12 GROUPS IN THE CHANELS GIVEN FOR EACH PATIENT
def group_inf(nomgroup, data_chan, upperchanels):
    group_date = []
    for e in nomgroup:
        if e in upperchanels:
            num = upperchanels.index(e)
            group_date.append(data_chan[num])
    return group_date


# DEF MEDIUM DELS GRUPS group_dat = [[1,2,3,..., 119345],[1,2,3,..., 119345],[1,2,3,..., 119345],[1,2,3,..., 119345]]
def mediumchanels(group_dat):
    segundo = len(group_dat[0])
    primero = len(group_dat)
    mean_data = []
    for i in range(0, segundo):
        suma1 = []
        for e in range(0, primero):
            suma1.append(group_dat[e][i])
        numpero = statistics.median(suma1)
        mean_data.append(numpero)
    return mean_data

# Obtencion de los grupos


# chanlesnames 62, groupp 12 i cada un amb els noms, data_chan all data 1 pacient EO o EC
def opteciogrups(chanlesnames, groupp, data_chan):
    group_date_finale = []  # Save all informatio chanels 12.
    group_date_final = []  # final copy
    upchan = upperchanel(chanlesnames)
    for i in groupp:  # x cada un dels 12 grups
        group_date = group_inf(i, data_chan, upchan)
        mean_group = mediumchanels(group_date)
        mean_array = np.array(mean_group)
        group_date_finale.append(mean_array)
    group_date_final = np.array(group_date_finale)
    return group_date_final


'''Execution'''

# IMPOT FUNTION THAT RETURNS 3 LISTS: A LIST OF ALL PAIENTS DATA, LIST NAMES PATIENS EO, LIST NAMES PATIENS EC:


def Import_Patients(basePATH):

    # DECLARE THE 12 CHANELS GROUPS WE WANT TO REDUCE THE 62 CHANELS
    Anterior_midline = ["FPZ", "AFZ", "FZ", "FCZ", "CZ"]
    Left_frontal = ["F7", "F5", "F3", "F1", "AF7", "AF3", "FP1"]
    Right_frontal = ["FP2", "AF4", "AF8", "F2", "F4", "F6", "F8"]
    Left_temporal = ["FT7", "T7", "TP7"]
    Left_central = ["FC5", "FC3", "FC1", "C5", "C3", "C1"]
    Left_parietal = ["CP5", "CP3", "CP1", "P7", "P5", "P3", "P1"]
    Left_Occipital = ["PO7", "PO3", "O1", "PO9"]
    Right_Occipital = ["PO4", "PO8", "O2", "PO10"]
    Right_parietal = ["P2", "P4", "P6", "P8", "CP2", "CP4", "CP6"]
    Right_temporal = ["FT8", "T8", "TP8"]
    Posterior_midline = ["CPZ", "PZ", "POZ", "OZ", "IZ"]
    Right_central = ["FC2", "FC4", "FC6", "C2", "C4", "C6"]

    # SAVE THIS 12 IN A MATRIX TO BE ITERATED LATER
    groups = [Anterior_midline, Left_frontal, Right_frontal, Left_temporal, Left_central, Left_parietal,
              Left_Occipital, Right_Occipital, Right_parietal, Right_temporal, Posterior_midline, Right_central]

    # READ FILES IN THE FOLDER AND SAVE THEM

    # Save the list of names separated EO vs EC
    listPATHEC = []
    listPATHEO = []
    # Save only the set files given that in the upload they automaticly search for the .fdt
    for root, dirs, files in os.walk(str(basePATH)):
        for filename in files:
            if '.set' in filename:
                if 'EC' in filename:
                    listPATHEC.append(filename)
                if 'EO' in filename:
                    listPATHEO.append(filename)
    # Save again but sorted by the name that way we know they are in order
    sorted_list_EO = sorted(listPATHEO)
    sorted_list_EC = sorted(listPATHEC)

    # ONCE WE HAVE THE NAMES IS TIME TO UPLOAD THE DATA
    # Variable were we will have all the patients with their EO and EC respective.
    EO_EC_Pacients = []

    for i in sorted_list_EO:
        # Load the doc
        x = mne.io.read_raw_eeglab(basePATH+i, preload=True, verbose=True)
        # GET DATA
        data = x._data
        # GET CHANELS
        chanles_names = x.ch_names
        # REDUCE CHANELS TO 12
        EO = opteciogrups(chanles_names, groups, data)

        # Select the same but with EC
        EC_name = [sorted_list_EC.index(e)
                   for e in sorted_list_EC if i[:-6] in e]
        EC_name_PATH = sorted_list_EC[EC_name[0]]
        # LOAD EC
        x2 = mne.io.read_raw_eeglab(
            basePATH+EC_name_PATH, preload=True, verbose=True)
        # GET DATA
        data2 = x2._data
        # GET CHANELS
        chanles_names2 = x2.ch_names
        # REDUCE CHANELS TO 12
        EC = opteciogrups(chanles_names2, groups, data2)

        # save
        EO_EC_P = []
        EO_EC_P.append(EO)
        EO_EC_P.append(EC)
        EO_EC_Pacients.append(EO_EC_P)

    return(sorted_list_EO, sorted_list_EC, EO_EC_Pacients)
