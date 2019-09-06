# coding: utf-8
#import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.pylab as pylab
from matplotlib.ticker import FormatStrFormatter
import os
import subprocess
import csv #to read the txt files
import sys
#import other self-defined functions
import __plotting_functions
import __fallspeed_relations
import __isa
import __postprocess_McSnow
#from IPython.core.debugger import Tracer ; Tracer()()
from matplotlib import rc

'''
this script compares the fitted data m/A-D with literature
'''


##general settings for the plotting
number_of_plots = 8

#optimize the appearance of the plot (figure size, fonts)
aspect_ratio=1./4.
legend_fontsize='medium'
#increase font sizes
params = {'legend.fontsize': legend_fontsize,
    'axes.labelsize': 'x-large', #size: Either a relative value of 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large' or an absolute font size, e.g., 12
    'axes.titlesize':'x-large',
    'xtick.labelsize':'x-large',
    'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
#define figure
figsize_height = 1./aspect_ratio*number_of_plots
fig, axes = plt.subplots(nrows=number_of_plots, ncols=1, figsize=(12.5,figsize_height))

#set up array of diameters for displaying the fits
low_diam_log= -4; high_diam_log=np.log10(4e-2)

#set up the diameter array
diam = np.logspace(low_diam_log,high_diam_log,500) #set up array of diameters
diam_log = np.log10(diam)
diam_center = (diam[:-1]+diam[1:])/2 #center of the diameter bins
diam_log_center = 10**((np.log10(diam[:-1])+np.log10(diam[1:]))/2) #logarithmic center of the diameter bins

####
#START: read modelled geometry (and vterm)
####
fit_dic = dict()
for key in ["particle_type","mass_allagg_coeff","area_allagg_coeff","mass_Nmono1_coeff","area_Nmono1_coeff"]:
    fit_dic[key] = []

with open("/home/mkarrer/Dokumente/plots/fit_results.txt","rb") as txtfile: #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in
    fit_param_reader = csv.reader(txtfile, delimiter='&', quoting=csv.QUOTE_NONE, lineterminator=os.linesep, escapechar=" ") #quoting avoids '' for formatted string; lineterminator avoids problems with system dependend lineending format https://unix.stackexchange.com/questions/309154/strings-are-missing-after-concatenating-two-or-more-variable-string-in-bash?answertab=active#tab-top
    for i_row,row_content in enumerate(fit_param_reader): #TODO: the row is not any more equal to the monomer number #read the row with N_mono_now (the currently considered monomer number)
        if row_content[0] in ["plate_prec","dendrite_prec","mixdendneedle_prec","needle_prec","column_prec","mixcolumndend_prec","bullet_prec","rosette_prec"] : #skip header and low precision data
            fit_dic["particle_type"] = np.append(fit_dic["particle_type"],row_content[0][:-5])
            fit_dic["mass_allagg_coeff"] = np.append(fit_dic["mass_allagg_coeff"],[float(row_content[-4]),float(row_content[-3])]) #python indices start with 0
            fit_dic["area_allagg_coeff"] = np.append(fit_dic["area_allagg_coeff"],[float(row_content[-2]),float(row_content[-1])])
            
            fit_dic["mass_Nmono1_coeff"] = np.append(fit_dic["mass_Nmono1_coeff"],[float(row_content[1]),float(row_content[4])]) #python indices start with 0
            fit_dic["area_Nmono1_coeff"] = np.append(fit_dic["area_Nmono1_coeff"],[float(row_content[7]),float(row_content[10])])

print fit_dic
#calculate the m/A-D array for each particle_type
for prop in ["mass","area"]:
    for i_particle_type,particle_type in enumerate(fit_dic["particle_type"]):
        fit_dic[prop + "_" + particle_type] = fit_dic[prop + "_allagg_coeff"][2*i_particle_type]*diam**fit_dic[prop + "_allagg_coeff"][2*i_particle_type+1]
        if prop=="area": #that means we got both mass and area already
            fit_dic["vterm_" + particle_type] = __fallspeed_relations.calc_vterm("bohm",fit_dic["mass" + "_" + particle_type],diam,fit_dic["area" + "_" + particle_type])
        fit_dic[prop + "_Nmono1_" + particle_type] = fit_dic[prop + "_Nmono1_coeff"][2*i_particle_type]*diam**fit_dic[prop + "_Nmono1_coeff"][2*i_particle_type+1]
        if prop=="area": #that means we got both mass and area already
            fit_dic["vterm_Nmono1_" + particle_type] = __fallspeed_relations.calc_vterm("bohm",fit_dic["mass_Nmono1_" + particle_type],diam,fit_dic["area_Nmono1_" + particle_type])
####
#END: read modelled geometry (and vterm)
####

#also show a sphere?
sphere_dic = dict()
rho_i = 917.6 #define ice density
sphere_dic["particle_type"] = "sphere"
sphere_dic["mass_allagg_coeff"] = [rho_i*np.pi/6.,3]
sphere_dic["area_allagg_coeff"] = [np.pi/4.,2]
for prop in ["mass","area"]:
    sphere_dic[prop] = sphere_dic[prop + "_allagg_coeff"][0]*diam**sphere_dic[prop + "_allagg_coeff"][1]
#calculate v-D
sphere_dic["vterm"] = __fallspeed_relations.calc_vterm("bohm",sphere_dic["mass"],diam,sphere_dic["area"])
####
#START: define some literature values
####

#######
#Mitchell 1996 #parameterization for m-D and A-D
#######
M96 = dict()
M96["particle_type"] = ["mixS3","sideplane","polycrys"]
M96["mass_coeff_mixS3"] = [0.0028*10.**(2.*2.1-3.),2.1,800e-6,4500e-6]
M96["area_coeff_mixS3"] = [0.2285*10.**(2.*1.88-4.),1.88] #area approximated as assemblage of planar polycrstals

M96["mass_coeff_sideplane"] = [0.0033*10.**(2.*2.2-3.),2.2,600e-6,4100e-6] 
M96["area_coeff_sideplane"] = [0.2285*10.**(2.*1.88-4.),1.88] #area approximated as assemblage of planar polycrstals

M96["mass_coeff_polycrys"] = [0.00739*10.**(2.*2.45-3.),2.45,20e-6,450e-6]
M96["area_coeff_polycrys"] = [0.2285*10.**(2.*1.88-4.),1.88]

#calculate array of geometric properties
for key in M96["particle_type"]:
    for prop in ["mass","area"]:
        M96[prop + '_' + key]= M96[prop + "_coeff_" + key][0]*diam**M96[prop + "_coeff_" + key][1]
    
#calculate v-D
for key in M96["particle_type"]:
    M96["vterm_" + key] = __fallspeed_relations.calc_vterm("bohm",M96["mass_" + key],diam,M96["area_" + key])

#####
#Locatelli & Hobbs 
#####
LH74 = dict()
LH74["particle_type"] = ["dendrite","mixed","sideplanes"]
LH74["mass_coeff_dendrite"] = [0.073*10.**(3*1.4-6.),1.4,2e-3,10e-3]
LH74["vterm_coeff_dendrite"] = [0.8*10.**(3*0.16),0.16,2e-3,10e-3]

LH74["mass_coeff_mixed"] = [0.037*10.**(3*1.9-6.),1.9,1e-3,3e-3]
LH74["vterm_coeff_mixed"] = [0.69*10.**(3*0.41),0.41,0.2e-3,3e-3]

LH74["mass_coeff_sideplanes"] = [0.04*10.**(3*1.4-6.),1.4,0.5e-3,4.0e-3]
LH74["vterm_coeff_sideplanes"] = [0.82*10.**(3*0.12),0.12,0.5e-3,4.0e-3]

#calculate mass and vterm array of geometric properties
for key in LH74["particle_type"]:
    for prop in ["mass","vterm"]:
        LH74[prop + '_' + key]= LH74[prop + "_coeff_" + key][0]*diam**LH74[prop + "_coeff_" + key][1]

#####
#Kajikawa 1989 
#####
K89 = dict()
K89["particle_type"] = ["dendriteP1e","stellarP2a","plateP2e"]

K89["mass_coeff_dendriteP1e"] = [4.82*10**(-4)*10.**(2*1.97-3.),1.97,0.08e-2,0.68e-2]
K89["vterm_mass_coeff_dendriteP1e"] = [102*10.**(2*0.11-2),0.11]

K89["mass_coeff_stellarP2a"] = [8.30*10**(-4)*10.**(2*2.09-3.),2.09,0.14e-2,0.70e-2]
K89["vterm_mass_coeff_stellarP2a"] = [94*10.**(2*0.089-2),0.089]

K89["mass_coeff_plateP2e"] = [3.96*10**(-4)*10.**(2*1.40-3.),1.40,0.08e-2,0.46e-2]
K89["vterm_mass_coeff_plateP2e"] = [68*10.**(2*0.053-2),0.053]

#calculate mass and vterm array of geometric properties
for key in K89["particle_type"]:
    for prop in ["mass"]:
        K89[prop + '_' + key]= K89[prop + "_coeff_" + key][0]*diam**K89[prop + "_coeff_" + key][1]
    for prop in ["vterm"]:
        K89[prop + '_' + key]= K89[prop + "_mass_coeff_" + key][0]*K89['mass_' + key]**K89[prop + "_mass_coeff_" + key][1]

#####
#Kajikawa 1973 individual crystals 
#####
K73 = dict()

K73["particle_type"] = ["plate","dendrite","thickplate"]
diamK73                 = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]) #[mm]
diamK73 = diamK73/1000 #[m]
K73["vterm_thickplate"] = np.array([10,  20, 30, 44, 57,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]) #[cm/s]
K73["vterm_plate"]      = np.array([np.nan,  10, 15, 20, 25,31,35,40,46,51,56,60,65,70,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]) #[cm/s]
K73["vterm_dendrite"]      = np.array([np.nan,  np.nan, np.nan,12,14,15,17,18,19,20,21,21.5,22,22.5,23,23.5,24,np.nan,  np.nan, np.nan]) #[cm/s]

for prop in K73["particle_type"]:
    K73["vterm_"+ prop] = K73["vterm_"+ prop]/100 #[m/s]


####
##Barthazy & Schefold
####
#constants for exponential fit v(D)=a(1-exp(-b D))
BS06 = dict()
BS06["particle_type"] = ["plate","modrimeddend","needle","modrimeddend","irreg"]

BS06["needle_expo"]=[1.10,2.86,0.298,2.682] #2 exponential fit parameter for aggregates of needles; range of the fit [mm]
BS06["plate_expo"]=[1.18,2.68,0.298,2.682] #2 exponential fit parameter for aggregates of needles; range of the fit
BS06["modrimeddend_expo"]=[1.29,1.58,0.447,3.725] #exponential fit parameter for aggregates of needles; range of the fit
BS06["irreg_expo"]=[1.57,1.95,0.149,4.768] #2 exponential fit parameter for aggregates of needles; range of the fit

for habit in BS06["particle_type"]:
    BS06["v_uncorr" + habit] = BS06[habit + "_expo"][0]*(1-np.exp(-BS06[habit + "_expo"][1]*(diam*1000)))
    rho_0m = __isa.isa(0)[2]; rho_MtRigi = __isa.isa(1604)[2]
    BS06["v_denscorr" + habit] = BS06["v_uncorr" + habit]*(rho_MtRigi/rho_0m)**0.54 #this is from Heymsfield (2007) as quoted in SZ10


####
##Brandes et. al 2008
####
#powerlaw fit for volume equivalent diameter
B08 = dict()
B08["particle_type"] = ["T-10","T-5","T-1"]

B08["vterm_coeff_T-10"] = [0.55*10.**(3*0.23),0.23,0.4e-3,14e-3] #size range: lower boundary: Brandes et. al 2007 Section 2 end; upper boundary: estimates from Fig. 3
B08["vterm_coeff_T-5"]=[0.67*10.**(3*0.25),0.25,0.4e-3,15e-3]
B08["vterm_coeff_T-1"]=[0.87*10.**(3*0.23),0.23,0.4e-3,20e-3]

for habit in B08["particle_type"]:
    for prop in ["vterm"]:
        B08[prop + '_' + habit]= B08[prop + "_coeff_" + habit][0]*diam**B08[prop + "_coeff_" + habit][1]
####
##Szyrmer and Zawadzki (2010)
####
#constant for the mean v(D)=au*Dmelted[cm]*bu
SZ10 = dict()

SZ10["particle_type"] = ["average"]

SZ10["mass_coeff_average"] = [0.0044*10.**(-3+2*2.0),2.0,1e-3,13e-3] #size range: lower boundary: Brandes et. al 2007 Section 2 end; upper boundary: estimates from Fig. 3

SZ10["vterm_coeff_average"] = [102*10.**(-2+2*0.18),0.23,1e-3,13e-3] #size range: lower boundary: Brandes et. al 2007 Section 2 end; upper boundary: estimates from Fig. 3
#SZ10["vterm_coeff_average"] = [2.337,0.18,1e-3,13e-3] #size range: lower boundary: Brandes et. al 2007 Section 2 end; upper boundary: estimates from Fig. 3

for habit in SZ10["particle_type"]:
    for prop in ["mass","vterm"]:
        SZ10[prop + '_' + habit]= SZ10[prop + "_coeff_" + habit][0]*diam**SZ10[prop + "_coeff_" + habit][1]
        
####
##von Lerber (2017)
####
#constant for the mean v(D)=au*Dmelted[cm]*bu
vL17 = dict()

vL17["particle_type"] = ["LWPsmall"] #,"LWPmedium","LWPlarge"]

vL17["mass_coeff_LWPsmall"] = [0.0046*10.**(-3+2*2.1),2.1,0.2e-3,13e-3] #size range: lower boundary: Section 1; upper boundary: taken from SZ10
vL17["mass_coeff_LWPmedium"] = [0.0053*10.**(-3+2*2.1),2.1,0.2e-3,13e-3] #size range: lower boundary: Section 1; upper boundary: taken from SZ10
vL17["mass_coeff_LWPlarge"] = [0.022*10.**(-3+2*2.3),2.3,0.2e-3,13e-3] #size range: lower boundary: Section 1; upper boundary: taken from SZ10

vL17["vterm_coeff_LWPsmall"] = [69*10.**(-2+2*0.20),0.20,0.2e-3,13e-3] #size range: lower boundary: Section 1; upper boundary: taken from SZ10
vL17["vterm_coeff_LWPmedium"] = [72*10.**(-2+2*0.24),0.24,0.2e-3,13e-3] #size range: lower boundary: Section 1; upper boundary: taken from SZ10
vL17["vterm_coeff_LWPlarge"] = [110*10.**(-2+2*0.33),0.33,0.2e-3,13e-3] #size range: lower boundary: Section 1; upper boundary: taken from SZ10

for habit in vL17["particle_type"]:
    for prop in ["mass","vterm"]:
        vL17[prop + '_' + habit]= vL17[prop + "_coeff_" + habit][0]*diam**vL17[prop + "_coeff_" + habit][1]
        
        
####
##von Lerber (2017) -with corrections from Annakaisa's mail from 4.8.
#### 
#constant for the mean v(D)=au*Dmelted[cm]*bu
vL17new = dict()

vL17new["particle_type"] = ["LWPsmall"] #,"LWPmedium","LWPlarge"]

vL17new["mass_coeff_LWPsmall"] = [0.0046*10.**(-3+2*2.1),2.1,0.2e-3,13e-3] #size range: lower boundary: Section 1; upper boundary: taken from SZ10
vL17new["mass_coeff_LWPmedium"] = [0.0053*10.**(-3+2*2.1),2.1,0.2e-3,13e-3] #size range: lower boundary: Section 1; upper boundary: taken from SZ10
vL17new["mass_coeff_LWPlarge"] = [0.022*10.**(-3+2*2.3),2.3,0.2e-3,13e-3] #size range: lower boundary: Section 1; upper boundary: taken from SZ10

vL17new["vterm_coeff_LWPsmall"] = [71.8*10.**(-2+2*0.20),0.20,0.2e-3,13e-3] #size range: lower boundary: Section 1; upper boundary: taken from SZ10
vL17new["vterm_coeff_LWPmedium"] = [75.5*10.**(-2+2*0.24),0.24,0.2e-3,13e-3] #size range: lower boundary: Section 1; upper boundary: taken from SZ10
vL17new["vterm_coeff_LWPlarge"] = [117.4*10.**(-2+2*0.33),0.33,0.2e-3,13e-3] #size range: lower boundary: Section 1; upper boundary: taken from SZ10

for habit in vL17new["particle_type"]:
    for prop in ["mass","vterm"]:
        vL17new[prop + '_' + habit]= vL17new[prop + "_coeff_" + habit][0]*diam**vL17new[prop + "_coeff_" + habit][1]
#other m-D 

other_mD = dict()

other_mD["particle_type"] = ["Cotton13","Baran11","Heym10"] #,"Cox88"]

#taken from Cotton et al (2013)
other_mD["mass_coeff_Cotton13"] = [0.0257,2.0,20e-6,500e-6] #in-situ aircraft
other_mD["mass_coeff_Baran11"] = [0.04,2.0,20e-6,500e-6] #forward modelling of reflectivity with prescribed exponent
other_mD["mass_coeff_Heym10"] = [0.176,2.2,20e-6,500e-6] #in-situ aircraft
other_mD["mass_coeff_Cox88"] = [0.069,2.0,20e-6,500e-6]


for habit in other_mD["particle_type"]:
    for prop in ["mass"]:
        other_mD[prop + '_' + habit]= other_mD[prop + "_coeff_" + habit][0]*diam**other_mD[prop + "_coeff_" + habit][1]
        
#calculate v-D with area from M96
for key in other_mD["particle_type"]:
    other_mD["vterm_" + key] = __fallspeed_relations.calc_vterm("bohm",other_mD["mass_" + key],diam,M96["area_mixS3"])
    
####
##Nettesheim et. al. (2019)
#### 
#constant for the mean v(D)=au*Dmelted[cm]*bu
Nett19 = dict()

Nett19["particle_type"] = ["plate","dendrite"] #,"LWPmedium","LWPlarge"]


Nett19["vterm_coeff_plate"] = [9.50e-1*10.**(6*0.13),0.13,-1.76,200e-6,5000e-6]
Nett19["vterm_coeff_dendrite"] = [1.94e-3*10.**(6*0.65),0.65,0.0,200e-6,5000e-6]

for habit in Nett19["particle_type"]:
    for prop in ["vterm"]:
        Nett19[prop + '_' + habit]= Nett19[prop + "_coeff_" + habit][0]*diam**Nett19[prop + "_coeff_" + habit][1]+Nett19[prop + "_coeff_" + habit][2]
        rho_ratio = 1000./900.#Nettesheim provided values for 900hPa #density scales with pressure 
        Nett19[prop + '_' + habit] = Nett19[prop + '_' + habit]*(rho_ratio)**0.54 #this is from Heymsfield (2007) as quoted in SZ10
####
#END: define some literature values
####
####
#START: read in Hyytiala data from Annakaisa
####

hyytiala_data = "/home/mkarrer/Dokumente/insitu_vterm/Hyytiala_MedianV_perc25_75.csv"
vD_hyyt_dic = dict() #dictionary which contains arrays of the observed variables from the Hyytiala site
vD_hyyt_dic["Dmax"] = []        #maximum dimension [mm]
vD_hyyt_dic["v_25perc"] = []    #25% percenticle of the velocity
vD_hyyt_dic["v_75perc"] = []    #75% percenticle of the velocity
vD_hyyt_dic["Dmax_wNAN"] = []   #maximum dimension, when including NaN velocity
vD_hyyt_dic["v"] = []           #median velocity
with open(hyytiala_data,"rb") as txtfile: #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in
    prop_reader = csv.reader(filter(lambda row: row[0]!='D',txtfile), delimiter=',', quoting=csv.QUOTE_NONNUMERIC, lineterminator=os.linesep) #row[0]!='#': skip the header; quoting avoids '' for formatted string; lineterminator avoids problems with system dependend lineending format https://unix.stackexchange.com/questions/309154/strings-are-missing-after-concatenating-two-or-more-variable-string-in-bash?answertab=active#tab-top
    #for i_row,row_content in enumerate(prop_reader):
    for row_content in prop_reader:
       
        #replace missing values and matlab "nan"'s with python "NaN"s
        for i_entry,entry in enumerate(row_content):
            if entry in ("","nan"):
                row_content[i_entry] = "NaN"

        vD_hyyt_dic["Dmax"] = np.append(vD_hyyt_dic["Dmax"],float(row_content[0]))
        vD_hyyt_dic["v_25perc"] = np.append(vD_hyyt_dic["v_25perc"],float(row_content[1]))
        vD_hyyt_dic["v_75perc"] = np.append(vD_hyyt_dic["v_75perc"],float(row_content[2]))
        vD_hyyt_dic["Dmax_wNAN"] = np.append(vD_hyyt_dic["Dmax_wNAN"],float(row_content[3]))
        vD_hyyt_dic["v"] = np.append(vD_hyyt_dic["v"],float(row_content[4]))

    
####
#END: read in Hyytiala data from Annakaisa
####

####
#START: read in CARE data from Annakaisa
####

hyytiala_data = "/home/mkarrer/Dokumente/insitu_vterm/CARE_2014_2017_V_D_median_25_75_larger05.csv"
vD_CARE_dic = dict() #dictionary which contains arrays of the observed variables from the Hyytiala site
vD_CARE_dic["Dmax"] = []        #maximum dimension [mm]
vD_CARE_dic["v_larger0"] = []           #median velocity
vD_CARE_dic["v_25perc_larger0"] = []    #25% percenticle of the velocity
vD_CARE_dic["v_75perc_larger0"] = []    #75% percenticle of the velocity
vD_CARE_dic["v_larger05"] = []           #median velocity
vD_CARE_dic["v_25perc_larger05"] = []    #25% percenticle of the velocity
vD_CARE_dic["v_75perc_larger05"] = []    #75% percenticle of the velocity
vD_CARE_dic["Dmax_wNAN"] = []   #maximum dimension, when including NaN velocity

with open(hyytiala_data,"rb") as txtfile: #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in
    prop_reader = csv.reader(filter(lambda row: row[0]!='D',txtfile), delimiter=';', quoting=csv.QUOTE_NONNUMERIC, lineterminator=os.linesep) #row[0]!='#': skip the header; quoting avoids '' for formatted string; lineterminator avoids problems with system dependend lineending format https://unix.stackexchange.com/questions/309154/strings-are-missing-after-concatenating-two-or-more-variable-string-in-bash?answertab=active#tab-top
    #for i_row,row_content in enumerate(prop_reader):
    for row_content in prop_reader:
       
        #replace missing values and matlab "nan"'s with python "NaN"s
        for i_entry,entry in enumerate(row_content):
            if entry in ("","nan"):
                row_content[i_entry] = "NaN"

        vD_CARE_dic["Dmax"] = np.append(vD_CARE_dic["Dmax"],float(row_content[0]))
        vD_CARE_dic["v_larger0"] = np.append(vD_CARE_dic["v_larger0"],float(row_content[1]))
        vD_CARE_dic["v_25perc_larger0"] = np.append(vD_CARE_dic["v_25perc_larger0"],float(row_content[2]))
        vD_CARE_dic["v_75perc_larger0"] = np.append(vD_CARE_dic["v_75perc_larger0"],float(row_content[3]))
        #vD_CARE_dic["Dmax"] = np.append(vD_CARE_dic["Dmax"],float(row_content[3]))
        vD_CARE_dic["v_larger05"] = np.append(vD_CARE_dic["v_larger05"],float(row_content[5]))
        vD_CARE_dic["v_25perc_larger05"] = np.append(vD_CARE_dic["v_25perc_larger05"],float(row_content[6]))
        vD_CARE_dic["v_75perc_larger05"] = np.append(vD_CARE_dic["v_75perc_larger05"],float(row_content[7]))
        #vD_CARE_dic["Dmax_wNAN"] = np.append(vD_CARE_dic["Dmax_wNAN"],float(row_content[3]))
    
####
#END: read in Hyytiala data from Annakaisa
####
###
#plot modelled and literature values
###

i_ax=0
for i_prop,(prop,prop_unit,prop_short) in enumerate(zip(["mass","mass_selected","area","vterm","vterm_frommass","vtermHyyt","vterm_selected","vterm_small"],["kg","kg","m2","ms-1","ms-1","ms-1","ms-1","ms-1"],["m","m","A","vterm","vterm","vterm","vterm","vterm"])):


    if prop in ["mass_selected","mass","area"]: #plot m/A-D

        #sphere
        if prop in ["mass_selected","mass"]:
            axes[i_ax].plot(diam,sphere_dic["mass"],linestyle="-",color="gray",label="sphere")
        elif prop=="area":
            axes[i_ax].plot(diam,sphere_dic[prop],linestyle="-",color="gray",label="sphere")
        #aggregate model fits
        for i_particle_type,particle_type in enumerate(fit_dic["particle_type"]):
            if not prop=="mass_selected":
                axes[i_ax].loglog(diam,fit_dic[prop + "_" + particle_type],linestyle=["-","--","--",":",":","--"][i_particle_type],marker=["","","o","","x","d"][i_particle_type],markevery=20,color=["g","g","g","g","g","limegreen"][i_particle_type],label="agg_model " + particle_type)
            else:
                axes[i_ax].loglog(diam,fit_dic["mass_" + particle_type],linestyle=["-","--","--",":",":","--"][i_particle_type],marker=["","","o","","x","d"][i_particle_type],markevery=20,color=["g","g","g","g","g","limegreen"][i_particle_type],label="agg_model " + particle_type)

        #M96
        for i_particle_type,particle_type in enumerate(M96["particle_type"]):
            diam_flagged = np.where(((diam>(M96["mass_coeff_" + particle_type][2])) & (diam<(M96["mass_coeff_" + particle_type][3]))),diam,np.nan)
            if prop in ["mass"]:
                axes[i_ax].loglog(diam_flagged,M96["mass_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="r",label="M96_" + particle_type)
            elif prop=="area":
                axes[i_ax].loglog(diam_flagged,M96[prop + "_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="r",label="M96_" + particle_type)
        if prop=="mass":
            #LH74
            for i_particle_type,particle_type in enumerate(LH74["particle_type"]):
                diam_flagged = np.where(((diam>(LH74["mass_coeff_" + particle_type][2])) & (diam<(LH74["mass_coeff_" + particle_type][3]))),diam,np.nan)
                axes[i_ax].loglog(diam_flagged,LH74[prop + "_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="orange",label="LH74_" + particle_type)
                
            #K89
            for i_particle_type,particle_type in enumerate(K89["particle_type"]):
                diam_flagged = np.where(((diam>(K89["mass_coeff_" + particle_type][2])) & (diam<(K89["mass_coeff_" + particle_type][3]))),diam,np.nan)
                axes[i_ax].loglog(diam_flagged,K89[prop + "_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="orangered",label="K89_" + particle_type)
            #vL17
            for i_particle_type,particle_type in enumerate(vL17["particle_type"]):
                diam_flagged = np.where(((diam>(vL17["mass_coeff_" + particle_type][2])) & (diam<(vL17["vterm_coeff_" + particle_type][3]))),diam,np.nan)
                axes[i_ax].plot(diam_flagged,vL17["mass_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="slateblue",label="vL17_" + particle_type)
            #SZ10
            for i_particle_type,particle_type in enumerate(SZ10["particle_type"]):
                diam_flagged = np.where(((diam>(SZ10["mass_coeff_" + particle_type][2])) & (diam<(SZ10["mass_coeff_" + particle_type][3]))),diam,np.nan)
                axes[i_ax].plot(diam_flagged,SZ10["mass_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="cyan",label="SZ10_" + particle_type)  
                
            for i_particle_type,particle_type in enumerate(other_mD["particle_type"]):
                diam_flagged = diam #np.where(((diam>(SZ10["mass_coeff_" + particle_type][2])) & (diam<(SZ10["mass_coeff_" + particle_type][3]))),diam,np.nan)
                axes[i_ax].plot(diam_flagged,other_mD["mass_" + particle_type],linestyle=["--","-.",":",":"][i_particle_type],color="magenta",label=particle_type)
        if prop=="mass_selected":
            #M96
            for i_particle_type,particle_type in enumerate(M96["particle_type"]):
                diam_flagged = np.where(((diam>(M96["mass_coeff_" + particle_type][2])) & (diam<(M96["mass_coeff_" + particle_type][3]))),diam,np.nan)
                axes[i_ax].loglog(diam_flagged,M96["mass_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="darkorange",label="M96_" + particle_type)
            #LH74
            for i_particle_type,particle_type in enumerate(LH74["particle_type"]):
                diam_flagged = np.where(((diam>(LH74["mass_coeff_" + particle_type][2])) & (diam<(LH74["mass_coeff_" + particle_type][3]))),diam,np.nan)
                axes[i_ax].loglog(diam_flagged,LH74["mass_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="gold",label="LH74_" + particle_type)
                
            for i_particle_type,particle_type in enumerate(["Cotton13"]):
                diam_flagged = np.where(((diam>(other_mD["mass_coeff_" + particle_type][2])) & (diam<(other_mD["mass_coeff_" + particle_type][3]))),diam,np.nan)
                axes[i_ax].plot(diam_flagged,other_mD["mass_" + particle_type],linestyle=["--","-.",":",":"][i_particle_type],color="magenta",label=particle_type)
    elif prop=="vterm":

        for scale in ["log"]: #["linear","log"]:
            #plot v-D
            
            #sphere
            axes[i_ax].plot(diam,sphere_dic["vterm"],linestyle="-",color="gray",label="sphere")
            
            #aggregate model fits
            for i_particle_type,particle_type in enumerate(fit_dic["particle_type"]):
                axes[i_ax].semilogx(diam,fit_dic[prop + "_" + particle_type],linestyle=["-","--","--",":",":","--"][i_particle_type],marker=["","","o","","x","d"][i_particle_type],markevery=20,color=["g","g","g","g","g","limegreen"][i_particle_type],label="agg_model " + particle_type)
            
            
            #LH74
            for i_particle_type,particle_type in enumerate(LH74["particle_type"]):
                diam_flagged = np.where(((diam>(LH74["vterm_coeff_" + particle_type][2])) & (diam<(LH74["vterm_coeff_" + particle_type][3]))),diam,np.nan)
                axes[i_ax].plot(diam_flagged,LH74["vterm_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="orange",label="LH74_" + particle_type)
                
            #K89
            for i_particle_type,particle_type in enumerate(K89["particle_type"]):
                diam_flagged = np.where(((diam>(K89["mass_coeff_" + particle_type][2])) & (diam<(K89["mass_coeff_" + particle_type][3]))),diam,np.nan)
                axes[i_ax].plot(diam_flagged,K89["vterm_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="orangered",label="K89_" + particle_type)
    elif prop=="vterm_selected":
        for scale in ["log"]: #["linear","log"]:
            #plot v-D
            
            #aggregate model fits
            for i_particle_type,particle_type in enumerate(fit_dic["particle_type"]):
                #if particle_type=="mixcolumndend": #dont show this sensitivity test here
                #    continue
                axes[i_ax].semilogx(diam,fit_dic["vterm_" + particle_type],linestyle=["-","--","--",":",":","--"][i_particle_type],marker=["","","o","","x","d"][i_particle_type],markevery=20,color=["g","g","g","g","g","limegreen"][i_particle_type],label="agg_model " + particle_type)
            #LH74
            for i_particle_type,particle_type in enumerate(LH74["particle_type"]):
                diam_flagged = np.where(((diam>(LH74["vterm_coeff_" + particle_type][2])) & (diam<(LH74["vterm_coeff_" + particle_type][3]))),diam,np.nan)
                axes[i_ax].plot(diam_flagged,LH74["vterm_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="orange",label="LH74_" + particle_type)
                
            
            #sphere
            axes[i_ax].plot(diam,sphere_dic["vterm"],linestyle="-",color="gray",label="sphere")
            
            #plot v-D quantiles from Hyytiala data
            axes[i_ax].plot(vD_hyyt_dic["Dmax"]*1e-3,vD_hyyt_dic["v"],linestyle="-",color="y",label="PIP Hyytiala(v>0.5m/s)")
            axes[i_ax].fill_between(vD_hyyt_dic["Dmax"]*1e-3,vD_hyyt_dic["v_25perc"],vD_hyyt_dic["v_75perc"],color="y",alpha=0.1)

            
            #plot v-D quantiles from CARE data
            axes[i_ax].plot(vD_CARE_dic["Dmax"]*1e-3,vD_CARE_dic["v_larger0"],linestyle="-",color="red",label="PIP CARE(v>0.0m/s)")
            axes[i_ax].fill_between(vD_CARE_dic["Dmax"]*1e-3,vD_CARE_dic["v_25perc_larger0"],vD_CARE_dic["v_75perc_larger0"],color="red",alpha=0.1)
            
            #plot v-D quantiles from CARE data
            axes[i_ax].plot(vD_CARE_dic["Dmax"]*1e-3,vD_CARE_dic["v_larger05"],linestyle="-",color="blue",label="PIP CARE(v>0.5m/s)")
            axes[i_ax].fill_between(vD_CARE_dic["Dmax"]*1e-3,vD_CARE_dic["v_25perc_larger05"],vD_CARE_dic["v_75perc_larger05"],color="blue",alpha=0.1)
            #i_ax+=1
            #i_ax+=1
        #i_ax+=-1
    elif prop=="vterm_small":
        for scale in ["log"]: #["linear","log"]:
            #plot v-D
            
            #aggregate model fits
            for i_particle_type,particle_type in enumerate(fit_dic["particle_type"]):
                if particle_type in ["mixcolumndend","mixdendneedle","needle"]: #dont show this sensitivity test here
                    continue
                axes[i_ax].semilogx(diam,fit_dic["vterm_" + particle_type],linestyle=["-","--",":",":",":"][i_particle_type],marker=["","","","","","d"][i_particle_type],markevery=20,color=["g","g","g","g","g"][i_particle_type],label="agg_model " + particle_type)
                axes[i_ax].semilogx(diam,fit_dic["vterm_Nmono1_" + particle_type],linestyle=["-","--","--",":",":",":"][i_particle_type],marker=["","","o","","","d"][i_particle_type],markevery=20,color=["b","b","b","b","b"][i_particle_type],label="agg_model mono " + particle_type)
            #K89
            #for i_particle_type,particle_type in enumerate(K89["particle_type"]):
            #    diam_flagged = np.where(((diam>(K89["mass_coeff_" + particle_type][2])) & (diam<(K89["mass_coeff_" + particle_type][3]))),diam,np.nan)
            #    axes[i_ax].loglog(diam_flagged,K89["vterm_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="limegreen",label="K89_" + particle_type)

            #K73
            for i_particle_type,particle_type in enumerate(K73["particle_type"]):
                axes[i_ax].loglog(diamK73,K73["vterm_" + particle_type],linestyle=["-","--",":",":",":"][i_particle_type],color="aqua",label="K73_" + particle_type)   

            #Nett19
            for i_particle_type,particle_type in enumerate(Nett19["particle_type"]):
                diam_flagged = np.where(((diam>(Nett19["vterm_coeff_" + particle_type][3])) & (diam<(Nett19["vterm_coeff_" + particle_type][4]))),diam,np.nan)
                axes[i_ax].loglog(diam_flagged,Nett19["vterm_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="magenta",label="Nett19_" + particle_type)                
            #sphere
            #axes[i_ax].plot(diam,sphere_dic["vterm"],linestyle="-",color="gray",label="sphere")
            
             
            #sphere
            axes[i_ax].plot(diam,sphere_dic["vterm"],linestyle="-",color="gray",label="sphere")


    elif prop=="vterm_frommass":
        for scale in  ["log"]: #["linear","log"]:
            ##sphere
            axes[i_ax].plot(diam,sphere_dic["vterm"],linestyle="-",color="gray",label="sphere")
            black cloth and it was dropped through the slit
(the outer cylinder was taken away and the shutter was opened beforehand).
            #aggregate model fits
            for i_particle_type,particle_type in enumerate(fit_dic["particle_type"]):
                axes[i_ax].semilogx(diam,fit_dic["vterm_" + particle_type],linestyle=["-","--","--",":",":","--"][i_particle_type],marker=["","","o","","x","d"][i_particle_type],markevery=20,color=["g","g","g","g","g","limegreen"][i_particle_type],label="agg_model " + particle_type)
            
            #M96
            for i_particle_type,particle_type in enumerate(["mixS3"]): #M96["particle_type"]):
                diam_flagged = diam # np.where(((diam>(M96["mass_coeff_" + particle_type][2])) & (diam<(M96["mass_coeff_" + particle_type][3]))),diam,np.nan)
                axes[i_ax].plot(diam_flagged,M96["vterm_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="r",label="M96_" + particle_type + "_bohm")
            
            #other: use several m-D fits and M96 area
            for i_particle_type,particle_type in enumerate(other_mD["particle_type"]):
                diam_flagged = diam #np.where(((diam>(other_mD["mass_coeff_" + particle_type][2])) & (diam<(other_mD["mass_coeff_" + particle_type][3]))),diam,np.nan)
                axes[i_ax].plot(diam_flagged,other_mD["vterm_" + particle_type],linestyle=["--","-.",":",":"][i_particle_type],color="r",label="M96_" + particle_type + "_bohm")
            #i_ax+=1
        #i_ax+=-1
    elif prop=="vtermHyyt":
        for scale in  ["log"]: #["linear","log"]:
            #sphere
            #axes[i_ax].plot(diam,sphere_dic["vterm"],linestyle="-",color="gray",label="sphere")
            
            #aggregate model fits
            #for i_particle_type,particle_type in enumerate(fit_dic["particle_type"]):
            #    axes[i_ax].semilogx(diam,fit_dic["vterm_" + particle_type],linestyle=["-","--","--",":",":","--"][i_particle_type],marker=["","","o","","x","d"][i_particle_type],markevery=20,color=["g","g","g","g","g","limegreen"][i_particle_type],label="agg_model " + particle_type)


            #BS06
            #for i_habit,habit in enumerate(BS06["particle_type"]): #"irreg"
            #    diam_flagged = np.where(((diam>(BS06[habit + "_expo"][2]*1e-3)) & (diam<(BS06[habit + "_expo"][3]*1e-3))),diam,np.nan)
            #    axes[i_ax].plot(diam_flagged,BS06["v_denscorr" + habit],linestyle=["-","--","-.",":",":"][i_habit],color="blue",label="BS06" + habit)    

            #B08
            #for i_particle_type,particle_type in enumerate(B08["particle_type"]):
            #    diam_flagged = np.where(((diam>(B08["vterm_coeff_" + particle_type][2])) & (diam<(B08["vterm_coeff_" + particle_type][3]))),diam,np.nan)
            #    axes[i_ax].plot(diam_flagged,B08["vterm_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="deepskyblue",label="B08_" + particle_type)

            #SZ10
            #for i_particle_type,particle_type in enumerate(SZ10["particle_type"]):
            #    diam_flagged = np.where(((diam>(SZ10["vterm_coeff_" + particle_type][2])) & (diam<(SZ10["vterm_coeff_" + particle_type][3]))),diam,np.nan)
            #    axes[i_ax].plot(diam_flagged,SZ10["vterm_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="cyan",label="SZ10_" + particle_type)
                
            #vL17
            #for i_particle_type,particle_type in enumerate(vL17["particle_type"]):
            #    diam_flagged = np.where(((diam>(vL17["vterm_coeff_" + particle_type][2])) & (diam<(vL17["vterm_coeff_" + particle_type][3]))),diam,np.nan)
            #    axes[i_ax].plot(diam_flagged,vL17["vterm_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="slateblue",label="vL17_" + particle_type)
            #vL17new
            #for i_particle_type,particle_type in enumerate(vL17new["particle_type"]):
            #    diam_flagged = np.where(((diam>(vL17new["vterm_coeff_" + particle_type][2])) & (diam<(vL17new["vterm_coeff_" + particle_type][3]))),diam,np.nan)
            #    axes[i_ax].plot(diam_flagged,vL17new["vterm_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="indigo",label="vL17new_" + particle_type)            
            #plot v-D quantiles from Hyytiala data
            #axes[i_ax].plot(vD_hyyt_dic["Dmax"]*1e-3,vD_hyyt_dic["v"],linestyle="-",color="y",label="Hyytiala")
            #axes[i_ax].fill_between(vD_hyyt_dic["Dmax"]*1e-3,vD_hyyt_dic["v_25perc"],vD_hyyt_dic["v_75perc"],color="y",alpha=0.3)
            #i_ax+=1
            #LH74
            for i_particle_type,particle_type in enumerate(LH74["particle_type"]):
                diam_flagged = np.where(((diam>(LH74["vterm_coeff_" + particle_type][2])) & (diam<(LH74["vterm_coeff_" + particle_type][3]))),diam,np.nan)
                axes[i_ax].plot(diam_flagged,LH74["vterm_" + particle_type],linestyle=["-","--","-.",":",":"][i_particle_type],color="y",label="LH74_" + particle_type)
                

            
            #plot v-D quantiles from Hyytiala data
            #axes[i_ax].plot(vD_hyyt_dic["Dmax"]*1e-3,vD_hyyt_dic["v"],linestyle="-",color="y",label="PIP Hyytiala(v>0.5m/s)")
            #axes[i_ax].fill_between(vD_hyyt_dic["Dmax"]*1e-3,vD_hyyt_dic["v_25perc"],vD_hyyt_dic["v_75perc"],color="y",alpha=0.1)

            
            #plot v-D quantiles from CARE data
            axes[i_ax].plot(vD_CARE_dic["Dmax"]*1e-3,vD_CARE_dic["v_larger0"],linestyle="-",color="g",label="PIP CARE(v>0.0m/s)")
            axes[i_ax].fill_between(vD_CARE_dic["Dmax"]*1e-3,vD_CARE_dic["v_25perc_larger0"],vD_CARE_dic["v_75perc_larger0"],color="g",alpha=0.1)
            
            #plot v-D quantiles from CARE data
            #axes[i_ax].plot(vD_CARE_dic["Dmax"]*1e-3,vD_CARE_dic["v_larger05"],linestyle="-",color="blue",label="PIP CARE(v>0.5m/s)")
            #axes[i_ax].fill_between(vD_CARE_dic["Dmax"]*1e-3,vD_CARE_dic["v_25perc_larger05"],vD_CARE_dic["v_75perc_larger05"],color="blue",alpha=0.1)

            #K73
            for i_particle_type,particle_type in enumerate(K73["particle_type"]):
                axes[i_ax].plot(diamK73,K73["vterm_" + particle_type],linestyle=["-","--",":",":",":"][i_particle_type],color="b",label="K73_" + particle_type)   
        #i_ax+=-1
    #do the labelling once in the end
    if prop in ["vterm","vterm_frommass","vtermHyyt","vterm_selected","vterm_small"]:
        for scale in ["log"]:
            #labelling
            axes[i_ax].set_xlabel("diameter D / m")
            axes[i_ax].set_ylabel(r"terminal velocity $v_{term}$ / m $s^{-1}$")
            axes[i_ax].set_xscale(scale)
            if prop in ["vterm_selected"]:
                axes[i_ax].legend(ncol=1,loc='center left', bbox_to_anchor=(1.0, 0.5))
            elif prop in ["vtermHyyt"]:
                axes[i_ax].legend(ncol=1)
            else:
                axes[i_ax].legend(ncol=2,loc='center left', bbox_to_anchor=(1.0, 0.5))

            axes[i_ax].grid(which="both")
            if scale=="linear":
                axes[i_ax].set_xlim([0,1e-2])
            else:
                axes[i_ax].set_xlim([1e-4,4e-2])


                
            axes[i_ax].set_ylim([0.0,2.0])
            if prop=="vterm_small": #overwrite limits
                axes[i_ax].set_xlim([1e-4,2e-3])
                axes[i_ax].set_ylim([0.0,1.0])
                axes[i_ax].set_yscale("linear")
            i_ax+=1
       
    else:

        axes[i_ax].grid(which="both")
        if prop in ["mass_selected"]:
            axes[i_ax].set_xlim([1e-4,1e-2])
        else:
            axes[i_ax].set_xlim([5e-4,1e-2])
        if prop in ["mass"]:
            axes[i_ax].set_ylim([1e-9,1e-5])
        elif prop in ["mass_selected"]:
            axes[i_ax].set_ylim([1e-11,1e-5])
        if prop=="area":
            axes[i_ax].set_ylim([1e-8,1e-4])
        #labelling
        axes[i_ax].set_xlabel("diameter D / m")
        if not prop=="mass_selected":
            axes[i_ax].set_ylabel((" / ").join((prop+ ' ' + prop_short,prop_unit)))
        else:
            axes[i_ax].set_ylabel((" / ").join(("mass " + prop_short,prop_unit)))
        #axes[i_ax].set_xscale(scale)
        if not prop in ["mass_selected"]:
            axes[i_ax].legend(ncol=2,loc='center left', bbox_to_anchor=(1.0, 0.5))
        else:
            axes[i_ax].legend(ncol=1,loc='center left', bbox_to_anchor=(1.0, 0.5))
        i_ax+=1


###########
###save the plot (and open it)
###########
plt.tight_layout()
dir_save = '/home/mkarrer/Dokumente/plots/geom_comp/'
out_filestring = "geom_modelvsobs"

plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
plt.savefig(dir_save + out_filestring + '.png', dpi=100)
print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])


ax_description_array=["m","m_selected","A","vterm_dirobs","vterm_indirobs","vterm_disdro","vterm_selected","vterm_small"]
for save_axis,ax_description in enumerate(ax_description_array):
    print "saving subfigure",save_axis,ax_description
    try: #clear label of axis before and recover current one
        axes[save_axis].set_xlabel("diameter D / m")
        axes[save_axis-1].set_xlabel("")
    except:
        pass
    extent = axes[save_axis].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    axes[save_axis-1].set_xticks([]) #remove previous xlabel so that is not visible
    axes[save_axis-1].set_xlabel("") #remove previous xlabel so that is not visible
    if ax_description=="vterm_disdro":
        fig.savefig('/home/mkarrer/Dokumente/plots/tmp.pdf',bbox_inches=extent.expanded(1.5, 1.4),dpi=400) #for the introduction it is nice to have the plot in the same size
    else:
        fig.savefig('/home/mkarrer/Dokumente/plots/tmp.pdf',bbox_inches=extent.expanded(2.3, 1.3),dpi=400) #(1.6, 1.4) are width and height expanded around center

    subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/tmp.pdf /home/mkarrer/Dokumente/plots/4paper/' + out_filestring + '_' + ax_description + '.pdf',shell=True)
    