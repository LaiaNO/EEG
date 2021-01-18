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
from numpy import asarray
from numpy import savetxt

def Import_Patients(basePATH):

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
    Chanels = []

    for i in sorted_list_EO:
        #Load the doc
        x=mne.io.read_raw_eeglab(basePATH+i, preload=True, verbose=True)
        #GET DATA
        EO = x._data
        #GET CHANELS
        chanles_names_EO = x.ch_names

        #Select the same but with EC
        EC_name = [sorted_list_EC.index(e) for e in sorted_list_EC if i[:-6] in e]
        EC_name_PATH = sorted_list_EC[EC_name[0]]
        #LOAD EC
        x2=mne.io.read_raw_eeglab(basePATH+EC_name_PATH, preload=True, verbose=True)
        #GET DATA
        EC = x2._data
        #GET CHANELS
        chanles_names_EC = x2.ch_names

        #save
        EO_EC_P = []
        EO_EC_P.append(EO)
        EO_EC_P.append(EC)
        EO_EC_Pacients.append(EO_EC_P)
        chan = []
        chan.append(chanles_names_EO)
        chan.append(chanles_names_EC)
        Chanels.append(chan)
    
    return(sorted_list_EO, sorted_list_EC, EO_EC_Pacients, Chanels)

#UPLOAD ALL THE FILES OF THE PATIENS IN A LIST LIKE THE FOLLOWING:
#Subject = EO_EC_Pacients[S] 
#EO = EO_EC_Pacients[S][0]
#EC = EO_EC_Pacients[S][1]
#Channel = EO_EC_Pacients[0][1][Num] 

#Subject = Chanels[S] 
#EO_Chanels_Names = Chanels[S][0]
#EC_Chanels_Names = Chanels[S][1]
#Channel_Name = Chanels[0][1][Num] 
basePATH = 'Files/Preprocessed/'
sorted_list_EO, sorted_list_EC, EO_EC_Pacients,Chanels = Import_Patients(basePATH)
print(len(sorted_list_EO))

