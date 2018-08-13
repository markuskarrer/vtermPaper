#import modulus
import numpy as np

# define functions to postprocess McSnow

def read_mass2frdat(experiment,filestring):
    '''
    read SP-properties from the mass2fr.dat file into the SP-dictionary
    INPUT experiment: descriptor (also file name) of the experiment
    '''
    directory = "/home/mkarrer/Dokumente/McSnow/MCSNOW/experiments/"
    #filestring = directory + experiment + '/mass2fr.dat'
    #load file from .dat
    SP_fullinfo = np.loadtxt(filestring)
    
    #create dictionary
    SP = dict()
    #names of variables in mass2fr
    varnames_mass2frdat = ["m_tot","Frim","height","d_rime","vt","xi",    "rhor","a_rime","mr_crit","diam",    "proj_A",   "mm"]
    #                       total   rime             ??    fall  multi-    rime     ?     crit mass  diameter  projected   monomer
    #                       mass    fraction         ??    speed plicity  density   ?  compl. infilling          area     multiplicity
    for i,key in enumerate(varnames_mass2frdat):
        SP[key] = SP_fullinfo[:,i]

    return SP

def seperate_by_height_and_diam(SP,nbins=100,nheights=51,model_top=500,diamrange=[-9,0]):
    '''
    count number of real particles (RP) in height-diameter bins
    INPUT:  SP:dictionary with properties of all SP of one timestep 
            nbins: number of bins in diameter array
            nheights: number of heights in height array
            model_top: top of model (also top of height area) [m]
            diamrange: specifies range of diameter array from 10**firstvalue to 10**secondvalue
    '''

    #create height vector
    zres = model_top/(nheights-1) #vertical resolution
    heightvec = np.linspace(0,model_top,nheights) #start with 0+zres and go n_heigts step up to model_top
    #define arrays with sp_diam
    d_bound_ds = np.logspace(diamrange[0],diamrange[1],nbins+1) #array from 1nm to 1m
    d_ds = d_bound_ds[:-1] + 0.5*np.diff(d_bound_ds) #diameter at center of bins
    #create dictionary for binned values and fill it with zeros
    binned_val = dict()
    binned_val["d_counts"] = np.zeros([nheights,nbins]) #number of RP at each h-D bin
    binned_val["d_counts_no_mult"] = np.zeros([nheights,nbins]) #number of SP at each h-D bin
    binned_val["av_Frim"] = np.zeros([nheights,nbins]) #mass averaged rime fraction for each h-D bin
    binned_val["av_rhor"] = np.zeros([nheights,nbins]) #mass averaged rime density for each h-D bin
    binned_val["av_mm"] = np.zeros([nheights,nbins]) #multiplicity weighted number
    #get the number of particles (SP*sp_multiplicity) at each bin (binned by height and sp_diameter)
    for i in range(0,nheights-1):
        for j in range(0,nbins):
                    condition_in_bin = np.logical_and(
                                    np.logical_and(d_bound_ds[j]<=SP["diam"],SP["diam"]<d_bound_ds[j+1]),
                                    np.logical_and(heightvec[i]<=SP["height"],SP["height"]<heightvec[i+1]),
                                    )
                    binned_val["d_counts"][i,j] = np.sum(np.where(condition_in_bin,SP["xi"],0))
                    binned_val["d_counts_no_mult"][i,j] = np.sum(np.where(condition_in_bin,1,0))
                    #get total number of RP per bin
                    multipl_bin = np.sum(np.where(condition_in_bin,SP["xi"],0))            
                    #get sum of qirim per bin
                    qirim_bin = np.sum(np.where(condition_in_bin,SP["m_tot"]*SP["Frim"]*SP["xi"],0))
                    #get sum of qitot
                    qitot_bin = np.sum(np.where(condition_in_bin,SP["m_tot"]*SP["xi"],0))
                    binned_val["av_Frim"][i,j] = qirim_bin/qitot_bin #calc. mass averaged rime fraction
                    #calc. rime mass averaged rime density
                    binned_val["av_rhor"][i,j] = np.sum(np.where(condition_in_bin,SP["m_tot"]*SP["Frim"]*SP["xi"]*SP["rhor"],0))/qitot_bin
                    #calc. multipl averaged
                    binned_val["av_mm"][i,j] = np.sum(np.where(condition_in_bin,SP["mm"]*SP["xi"],0))/multipl_bin

    binned_val["RPpSP"] = binned_val["d_counts"]/binned_val["d_counts_no_mult"]

    return binned_val,heightvec,d_bound_ds,d_ds,zres

def calculate_Vbox(experiment,zres):
    '''
    calculate the volume of the self (considering h) and the McSnow defined (considering A) box volume
    INPUT:  experiment: experiment-descriptor which includes information on runscripts
            delta_z: height steps in self-defined z-vector
    '''
    import re #for searching and selecting in strings
    #read box area from string
    box_area = float(re.search('ba(.*)', experiment).group(1)) *1./100. #this line gets the values after ba and before /mass2fr ; 1/100m^2 is the specification of the unit in the McSnow runscripts
    Vbox=box_area*zres
    return Vbox
def conv_num2numpm3(var,Vbox):
    #convert number in height bin [#] to numper per volume [#/m3]
    return var/Vbox
    