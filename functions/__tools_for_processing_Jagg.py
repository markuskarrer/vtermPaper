'''
some outsourced functions for the processing of Jussis-aggregates
'''
import numpy as np
from numpy import pi, r_
import matplotlib.pyplot as plt
from scipy import optimize
import sys

def init_particle_dict():
    '''
    initialize the particle dictionary
    '''
    particle_dic = dict() #dictionary which contains information about the type of the pristine crystals (dentrites, needles, plates, ...) the underlying distribution of the pristine crystals and the properties of the aggregates (m,A,D,...)
    particle_dic["particle_type"] = [] #type of the pristine crystals
    particle_dic["N_monomer"] = [] #number of monomers
    particle_dic["mass"] = [] #particle mass
    particle_dic["area"] = [] #particle area
    particle_dic["diam"] = [] #particle diameter
    #particle_dic["HWJussi"] = [] #particle area
    #particle_dic["KCJussi"] = [] #particle diameter
    
    return particle_dic

def fit_data(xdata,ydata,func='powerlaw',method='leastsq'):
    '''source: https://scipy-cookbook.readthedocs.io/items/FittingData.html
    fit any 2-dimensional data with different methods and target functions
    xdata: ordinate points
    ydata: coordinate points
    func: select the function to which the data should be fitted
    '''


    if func=='powerlaw':
        # Power-law fitting is best done by first converting
        # to a linear equation and then fitting to a straight line.
        # Note that the `logyerr` term here is ignoring a constant prefactor.
        #
        #  y = a * x^b
        #  log(y) = log(a) + b*log(x)
        #

        
        x = np.log10(xdata)
        y = np.log10(ydata)
        yerr = 2.0*  np.ones_like(ydata) #yerr / ydata  #as long as it is constant this doesnt matter

        # define our (line) fitting function
        fitfunc = lambda p, x: p[0] + p[1] * x
        errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

        #guess fit results first
        pinit = [1.0, 1.0]
    elif func=='Atlas':
        #Atlas-type fittong
        # y=A-B*exp(-gam*x)
        
        #stay linear for fitting
        x = (xdata)
        y = (ydata)
        
        yerr = 2.0*  np.ones_like(ydata) #yerr / ydata  #as long as it is constant this doesnt matter

        #define fit function
        fitfunc = lambda p, x: p[0] - p[1] * np.exp(-p[2]*x)
        errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
        #guess fit results first
        pinit = [1.2, 1.0,1e3]
        
    else: 
        print "function: " + function + " not implemented"
        sys.exit(1)
        
    if method=='leastsq':
        out = optimize.leastsq(errfunc, pinit,
                        args=(x, y, yerr), full_output=1)
    else: 
        print "method: " + method + " not implemented"
        sys.exit(1)
        
    #TODO: how is this done
    pfinal = out[0]
    covar = out[1]
    #print covar
    
    if func=='powerlaw': #convert the fitted variables back 
        #pfinal[ = pfinal[1]
        pfinal[0] = 10.0**pfinal[0]

    
    return pfinal,covar # amp,index


def fitting_wrapper(ax,particle_dic,prop,N_mono_list,usecolors,fit_dic,function="powerlaw"):
    '''
    select a subset of the particle in the particle_dic and
    fit functions (e.g. power law, Atlas type) to a certain property
    INPUT:
        - ax: axes on which the fitted function is added
        - particle_dic: dictionary containing the information from the aggregate property file
        - prop: name of the property which is fitted against the diameter (e.g. mass, area, fall speed)
        - N_mono_list: a list of the considered monomer numbers
        - usecolors: list of colors (which defines the colors for each monomer number)
        - fit_dic: dictionary which contains the fit-coefficients and the fitted values (e.g array of masses and corresponding diameter)
        - function: the functional relationship to which the property is fitted (e.g. powerlaw (default), Atlastype)
    OUTPUT:
        - ax: the axis on which the fits are added as lines
    
    '''
    ####
    #fit a power-law for each quantitity and Atlas-type function for "vterm_" 
    ####


    #if prop.startswith("vterm_"): #indirect fitting via the m-D and A-D coefficients
    #    fitresult_indir_powerlaw_string  = 'fit a*D**b (indirect)     a                 b\n' #string which contains the fit results
    
    #define some overall properties of the fit-dictionary and initialize it
    #set up a diamter array (for fitting indirectly (over the m-D,A-D coefficients) and plotting)
    diam = np.logspace(-4,-1,100)
    fit_dic["diam"] = diam
    fit_dic["N_mono_list"] = N_mono_list

    for i_N_mono,N_mono_now in enumerate(N_mono_list): #loop over all monomer numbers
        if not any(particle_dic["N_monomer"]==N_mono_now): #if there are no particles with this monomer number just skip this loop
            break
        
        # make a fit for the current number of monomers (power-law; always directly to the points coming from the aggregate model)
        [fitresult,covar] = fit_data(particle_dic["diam"][particle_dic["N_monomer"]==N_mono_now],particle_dic[prop][particle_dic["N_monomer"]==N_mono_now],func=function)
            
        
        
        #if prop.startswith("vterm_"): #fit additionally an Atlas-type function
        #    [fitresultAtlas,covar] = __tools_for_processing_Jagg.fit_data(particle_dic["diam"][particle_dic["N_monomer"]==N_mono_now],particle_dic[prop][particle_dic["N_monomer"]==N_mono_now],func='Atlas')
        
        #if prop=="mass":
        fit_dic[prop + "_coeff_Nmono_" + str(N_mono_now)] = fitresult #fitresult_arr_dir_powerlaw_string[0,i_N_mono] = fitresult #copy the A-D coefficients in the array
        fit_dic[prop + "_Nmono_" + str(N_mono_now)] = fit_dic[prop + "_coeff_Nmono_" + str(N_mono_now)][0]*fit_dic["diam"]**fit_dic[prop + "_coeff_Nmono_" + str(N_mono_now)][1] #calculate arrays of masses and areas based on the m-D fits
        #elif prop=="area":
        #    fit_dic[prop"coeff"] = fitresult#fitresult_arr_dir_powerlaw_string[1,i_N_mono] = fitresult #copy the A-D coefficients in the array
        
        
        #elif prop.startswith("vterm_"): #calculate also the fall speeds for an equidistant array corresponding to the m-D, A-D fit
        #    #calculate the masses and areas based on the fits of these properties and the previously defined diameter array
        #    mass_fitted = fitresult #fitresult_arr_dir_powerlaw_string[0,i_N_mono,0]*diam**fitresult_arr_dir_powerlaw_string[0,i_N_mono,1]
        #    area_fitted = fitresult #fitresult_arr_dir_powerlaw_string[1,i_N_mono,0]*diam**fitresult_arr_dir_powerlaw_string[1,i_N_mono,1]
            
        #calculate the fall speed based on the A-D and m-D fits
        #velocity_model = prop[6:] #get the name of the fall speed model from the property name (prop) which starts with "vterm"
        #particle_dic["vterm_" + velocity_model + "_fitted_via_mD_AD"] = __fallspeed_relations.calc_vterm(velocity_model,mass_fitted,diam,area_fitted)
    
        # make a fit for the current number of monomers (power-law; always directly to the points coming from the aggregate model)
        #[fitresult_indirect,covar] = __tools_for_processing_Jagg.fit_data(diam,particle_dic["vterm_" + velocity_model + "_fitted_via_mD_AD"],func='powerlaw')
        
        
        
        # output and saving of the fit results
        #print "fit coefficients (powerlaw): ",N_mono_now, " : ", fitresult
        #if prop.startswith("vterm_"): #fit additionally an Atlas-type function
        #    print "fit coefficients (Atlas): ",N_mono_now, " : ", fitresultAtlas

        #plot the power-law 
        ax.plot(diam,fitresult[0]*diam**fitresult[1],color=usecolors[i_N_mono],linewidth=0.5)



        #if prop.startswith("vterm_"): #add a line for the indirect fitting via the m-D and A-D coefficients
        #    pass #ax.plot(diam,fitresult_indirect[0]*diam**fitresult_indirect[1],color=usecolors[i_N_mono],linewidth=0.5,linestyle='--')
        #TODO: plot the _fitted_via_mD_AD data
        #if prop.startswith("vterm_"):
        #    #plot the scatter plot for the indirect data # use the colormap and the bounds which where defined before
        #    indirect_vD = ax.scatter(diam,particle_dic["vterm_" + velocity_model + "_fitted_via_mD_AD"],s=5,linewidth=1,facecolor=usecolors[i_N_mono],rasterized=True,marker='x') #,norm=colors.LogNorm(),cmap="gist_ncar")
    
    #now make a fit for all particle in the N_mono_list
        
        # make a fit for all aggregates in N_mono_list except N_mono=1 (power-law; always directly to the points coming from the aggregate model)
        N_mono_in_list = [False]*particle_dic["N_monomer"].shape[0]
        for N_mono_now in N_mono_list[1:]:
            N_mono_matching_now = (N_mono_now==particle_dic["N_monomer"])
            N_mono_in_list = np.logical_or(N_mono_matching_now,N_mono_in_list)

        [fitresult,covar] = fit_data(particle_dic["diam"][N_mono_in_list],particle_dic[prop][N_mono_in_list],func='powerlaw')
        fit_dic[prop + "_coeff_Nmono_allagg"] = fitresult
        fit_dic[prop + "_Nmono_allagg"] = fit_dic[prop + "_coeff_Nmono_allagg"][0]*fit_dic["diam"]**fit_dic[prop + "_coeff_Nmono_allagg"][1] #calculate arrays of masses and areas based on the m-D fits
        #plot the power-law 
        ax.plot(diam,fitresult[0]*diam**fitresult[1],color='black',linewidth=0.5)
        
        
    return ax, particle_dic #particle_dic is returned, because vterm_..._calculated via mD_AD is added
    
    
def add_fitresult_text(ax,fit_dic,prop,N_mono_list,function="powerlaw"):
    '''
    adds the coefficients of the fit result in the form of a text
    INPUT:
        ax: the axis on which the text is added (this is also returned)
        fit_dic: a dictionary which contains the coefficients of the fits
    '''
    
    #initialize a string with the result of the powerlaw and  Atlas-type fits
    if function=="powerlaw":
        fitresult_string = 'fit a*D**b         a                 b\n' #string which contains the fit results
    elif function=="Atlas":
        fitresult_string =    'fit A-B*exp(-C*D) A                B            C\n' #string which contains the fit results
    for N_mono_now in N_mono_list:
        fitresult_string += "Nmono=" + str(N_mono_now) + ': ' + str(fit_dic[prop  + "_coeff_Nmono_" + str(N_mono_now)]) + '\n' #write the current fit result into the general string
    fitresult_string += "Nmono>1: " + str(fit_dic[prop  + "_coeff_Nmono_allagg"]) + '\n'  #add also the result for all monomer numbers
    
    
    #if prop.startswith("vterm_"):  #indirect fitting via the m-D and A-D coefficients
    #    #fitresult_indir_powerlaw_string += "Nmono=" + str(N_mono_now) + ': ' + str(fitresult_indirect) + '\n' #write the current fit result into the general string
    #    fitresult_dir_Atlas_string += "Nmono=" + str(N_mono_now) + ': ' + str(fitresultAtlas) + '\n' #write the current fit result into the general string
    #    #plot the Atlas-type fit
    #    #ax.plot(diam,fitresultAtlas[0]-fitresultAtlas[1]*np.exp(-fitresultAtlas[2]*diam),color=usecolors[i_N_mono],linewidth=0.5,linestyle='-.')
    
    #add the text of the fitresult to the plot
    ax.text(0.0,1.0,fitresult_string,horizontalalignment='left',verticalalignment='top',transform = ax.transAxes,fontsize=4) #transform: give number in relative axes units instead of data coords
    #if prop.startswith("vterm_"):
    #    ax.text(0.0,0.6,fitresult_string,horizontalalignment='left',verticalalignment='top',transform = ax.transAxes,fontsize=4) #transform: give number in relative axes units instead of data coords

    #if prop.startswith("vterm_"):  #indirect fitting via the m-D and A-D coefficients
    #    pass #ax.text(0.0,0.6,fitresult_indir_powerlaw_string,horizontalalignment='left',verticalalignment='top',transform = ax.transAxes,fontsize=4) #transform: give number in relative axes units instead of data coords
    
    return ax