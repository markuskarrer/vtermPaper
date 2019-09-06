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
this script plots Annakaisa's measurement and adds older literature data on top of the plot
'''


##general settings for the plotting
number_of_plots = 3

#optimize the appearance of the plot (figure size, fonts)
[fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots,legend_fontsize='medium')


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

i_ax=1
for scale in ["linear","log"]:
    #plot v-D
    axes[i_ax].plot(vD_hyyt_dic["Dmax"]*1e-3,vD_hyyt_dic["v"],linestyle="-",color="g",label="Hyytiala")
    axes[i_ax].fill_between(vD_hyyt_dic["Dmax"]*1e-3,vD_hyyt_dic["v_25perc"],vD_hyyt_dic["v_75perc"],color="g",alpha=0.3)
    i_ax+=1
    
####
#END: read in Hyytiala data from Annakaisa
####


#set up array of diameters for displaying the fits
low_diam_log= -4; high_diam_log=np.log10(4e-2)

#set up the diameter array
diam = np.logspace(low_diam_log,high_diam_log,500) #set up array of diameters
diam_log = np.log10(diam)
diam_center = (diam[:-1]+diam[1:])/2 #center of the diameter bins
diam_log_center = 10**((np.log10(diam[:-1])+np.log10(diam[1:]))/2) #logarithmic center of the diameter bins




####
#START: overlay some literature values
####

'''
####
##Szyrmer and Zawadzki (2010)
####
#constant for the mean v(D)=au*Dmelted[cm]*bu
au=102 #[cm^0.82 s-1]
bu=0.18
#constant for mean m=am*Dmax*bm
am=0.0044 #[g cm-2]
bm=2.
rho_i = 917.6 #[kg m-3] #define ice density
#convert the maximum dimension array to the melted diameter in cm
diam_melted=((6.*1e-3*am*(diam*100)**bm)/(np.pi*rho_i))**(1./3.) #[m] #melted diameter corresponding to the above defined diam-array #ATTENTION: strange units: (6.*1e-3g->kg*am*(diam*100.cm->m)**bm)
diam_melted_flag = np.where((0.01/100<diam_melted) & (diam_melted<0.3/100)) #,diam_melted,np.nan)
v_SZ10 = au*(diam_melted*100)**bu/100. #[m/s] #au*diam_melted**bu/100. #[m/s]

#whats wrong here??
#axes[0].plot(diam_melted*100,v_SZ10,linestyle="--",color="b")
axes[0].plot(diam_melted[diam_melted_flag]*100,v_SZ10[diam_melted_flag],linestyle="-",color="b")
axes[0].set_xlabel("D_melted [cm]")

i_ax=1
for scale in ["linear","log"]:
    #plot v-D
    axes[i_ax].plot(diam,v_SZ10,linestyle="--",color="b")
    axes[i_ax].plot(diam[diam_melted_flag],v_SZ10[diam_melted_flag],linestyle="-",color="b",label="SZ10")
    i_ax+=1
'''



####
##Barthazy & Schefold
####
#constants for exponential fit v(D)=a(1-exp(-b D))
Barth = dict()

Barth["needles_expo"]=[1.10,2.86,0.298,2.682] #2 exponential fit parameter for aggregates of needles; range of the fit [mm]
Barth["plates_expo"]=[1.18,2.68,0.298,2.682] #2 exponential fit parameter for aggregates of needles; range of the fit
Barth["irreg_expo"]=[1.57,1.95,0.149,4.768] #2 exponential fit parameter for aggregates of needles; range of the fit
Barth["modrimeddend_expo"]=[1.29,1.58,0.447,3.725] #exponential fit parameter for aggregates of needles; range of the fit

for habit in ["plates","needles","irreg","modrimeddend"]:
    Barth["v_uncorr" + habit] = Barth[habit + "_expo"][0]*(1-np.exp(-Barth[habit + "_expo"][1]*(diam*1000)))
    rho_0m = __isa.isa(0)[2]; rho_MtRigi = __isa.isa(1604)[2]
    Barth["v_denscorr" + habit] = Barth["v_uncorr" + habit]*(rho_MtRigi/rho_0m)**0.54 #this is from Heymsfield (2007) as quoted in SZ10

i_ax=1
for scale in ["linear","log"]:
    #plot v-D
    for i_habit,habit in enumerate(["plates","modrimeddend","needles"]): #"irreg"
        #axes[i_ax].plot(diam,Barth["v_denscorr" + habit],linestyle=["--","-.",":"][i_habit],color="g",linewidth=)
        diam_flagged = np.where(((diam>(Barth[habit + "_expo"][2]*1e-3)) & (diam<(Barth[habit + "_expo"][3]*1e-3))),diam,np.nan)
        axes[i_ax].plot(diam_flagged,Barth["v_denscorr" + habit],linestyle=["--","-.",":"][i_habit],color="orange",label="BS05" + habit)

    i_ax+=1

'''
####
##von Lerber et. al. (2017)
####
vLerb17 = dict()
vLerb17["LWPsmall"] = [69./100,0.2] #[cm^(1-bv)s-1 1]

vLerb17["v"] = vLerb17["LWPsmall"][0]*(diam*100)**vLerb17["LWPsmall"][1]
#from IPython.core.debugger import Tracer ; Tracer()()

i_ax=1
for scale in ["linear","log"]:
    #plot v-D
    axes[i_ax].plot(diam,vLerb17["v"],linestyle="-",color="orange",label="vL17")

    i_ax+=1
'''

###
#Locatelli& Hobbs (1974)
###
LH74 = dict()
LH74["unr_side_planes"] = [0.81,0.99,0.4e-3,1.2e-3]
LH74["unr_agg_side_planes"] = [0.82,0.12,0.5e-3,4e-3] 
LH74["unr_agg_dend"] = [0.8,0.16,2.0e-3,12e-3] 
LH74["unr_agg_mix"] = [0.69,0.41,0.2e-3,3e-3]

habits = ["unr_side_planes","unr_agg_side_planes","unr_agg_dend","unr_agg_mix"]
for key in habits: #calculate vterm for each habit
    LH74["v_" + key] = LH74[key][0]*(diam*1000)**LH74[key][1]

i_ax=1
for scale in ["linear","log"]:
    #plot v-D
    for i_key,key in enumerate(habits):
        diam_flagged = np.where(((diam>(LH74[key][2])) & (diam<(LH74[key][3]))),diam,np.nan)
        axes[i_ax].plot(diam_flagged,LH74["v_" + key],linestyle=np.array(["-","--","-.",":"])[i_key],color="magenta",label="LH74_" + key)
    i_ax+=1

i_ax=1
    
#######
#Mitchell 1996 #parameterization for m-D and A-D
#######
M96 = dict()
M96["mixedagg_mass_coeff"] = [0.0028*10.**(2.*2.1-3.),2.1,800e-6,4500e-6]
M96["mixedagg_area_coeff"] = [0.2285*10.**(2.*1.88-4.),1.88]

#calculate area of geometric properties
for prop in ["mass","area"]:
    M96["mixedagg_" + prop]=M96["mixedagg_" + prop + "_coeff"][0]*diam**M96["mixedagg_" + prop + "_coeff"][1]
    
#calculate v-D
M96["vterm_mixedagg"] = __fallspeed_relations.calc_vterm("bohm",M96["mixedagg_mass"],diam,M96["mixedagg_area"])

'''
#mth,unr_alf,unr_bet,rhoi,rhol,Dth,unr_sig,unr_gam,sph_sig,sph_gam =  __postprocess_McSnow.return_parameter_mD_AD_rel("1d__xi")

M96["mixedagg_mass"]=unr_alf*diam**unr_bet
M96["mixedagg_area"]=unr_sig*diam**unr_gam
print unr_alf,unr_bet,unr_sig,unr_gam
print M96["mixedagg_mass_coeff"],M96["mixedagg_area_coeff"]
M96["vterm_mixedagg"] =  __fallspeed_relations.calc_vterm("bohm",M96["mixedagg_mass"],diam,M96["mixedagg_area"])
'''
i_ax=1
for scale in ["linear","log"]:
    #plot v-D
    diam_flagged = np.where(((diam>(M96["mixedagg_mass_coeff"][2])) & (diam<(M96["mixedagg_mass_coeff"][3]))),diam,np.nan)
    axes[i_ax].plot(diam_flagged,M96["vterm_mixedagg"],linestyle=":",color="k",label="M96_bohm_mixedagg")
    i_ax+=1



i_ax=1
#do the labelling once in the end
for scale in ["linear","log"]:
    #labelling
    axes[i_ax].set_xlabel("diameter D / m")
    axes[i_ax].set_ylabel(r"terminal vel. $v_{term}$ / m $s^{-1}$")
    axes[i_ax].set_xscale(scale)
    i_ax+=1
####
#END: overlay some literature values
####


###########
###save the plot (and open it)
###########
dir_save = '/home/mkarrer/Dokumente/plots/insitu_vterm/'
out_filestring = "vterm_insitu"


for ax in axes:
    ax.grid(which="both")
    ax.set_xlim([1e-4,4e-2])
    ax.set_ylim([0.0,2.5])
    plt.legend(ncol=3)
plt.tight_layout()

plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
plt.savefig(dir_save + out_filestring + '.png', dpi=100)
print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])

# Save just the portion _inside_ the .. axis's boundaries
save_axis=2
for i_ax in range(0,len(axes)):#clear all other axes
    if i_ax==save_axis:
        continue
    axes[i_ax].set_xlabel("")
    axes[i_ax].axes.get_xaxis().set_ticklabels([])
extent = axes[save_axis].get_window_extent().transformed(fig.dpi_scale_trans.inverted())

fig.savefig('/home/mkarrer/Dokumente/plots/tmp.pdf',bbox_inches=extent.expanded(1.5, 1.4),dpi=400)

subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/tmp.pdf' + ' ' + '/home/mkarrer/Dokumente/plots/4paper/' + out_filestring + '.pdf',shell=True)

