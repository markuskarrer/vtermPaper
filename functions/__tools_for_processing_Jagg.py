'''
some outsourced functions for the processing of Jussis-aggregates
'''
import numpy as np
from numpy import pi, r_
import matplotlib.pyplot as plt
from scipy import optimize
import sys
#from IPython.core.debugger import Tracer ; Tracer()()

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
    particle_dic["sizeparam_index"] = [] #size parameter
    #particle_dic["HWJussi"] = [] #particle area
    #particle_dic["KCJussi"] = [] #particle diameter
    
    return particle_dic

def fit_data(xdata,ydata,func='powerlaw',method='leastsq',weight='None'):
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
        yerr = 2.0*  np.ones_like(ydata)  #*  1./weight #* 1./ydata #yerr / ydata  #as long as it is constant this doesnt matter

        # define our (line) fitting function
        fitfunc = lambda p, x: p[0] + p[1] * x
        errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

        #guess fit results first
        pinit = [np.log10(0.01), 1.5]
    elif func=='mod_powerlaw':
        # Power-law fitting is best done by first converting
        # to a linear equation and then fitting to a straight line.
        # Note that the `logyerr` term here is ignoring a constant prefactor.
        #
        #  in this modified version of the power-law we let both coefficients vary linearly with the monomer number
        #  y = a *(1+a_mod*N_mono)*x^(b*(1+b_mod*N_mono)
        #  log(y) = log(a*((1+a_mod*N_mono)) + b*(1+b_mod*N_mono)*log(x)
        #
        
        x = np.log10(xdata[0])
        x2 = xdata[1]

        y = np.log10(ydata)
        yerr = 2.0*  np.ones_like(ydata)  #*  1./weight #* 1./ydata #yerr / ydata  #as long as it is constant this doesnt matter

        # define our (line) fitting function #inspired by: https://stackoverflow.com/questions/29374296/two-dimensional-fit-with-python
        fitfunc = lambda p, x, x2: p[0]*(1+p[1]*(x2-1)) + p[2]*(1+p[3]/(x2))* x
        errfunc = lambda p, x, x2, y, err: (y - fitfunc(p, x, x2)) / err

        #guess fit results first
        pinit = [np.log10(0.01), 0.2,2.5,-0.2]


    elif func=='Atlas':
        #Atlas-type fitting
        # y=A-B*exp(-gam*x)
        
        #stay linear for fitting
        x = (xdata)
        y = (ydata)
        
        #from IPython.core.debugger import Tracer ; Tracer()()
        if weight=='None': #dont apply any weighting
            yerr = 2.0 * np.ones_like(ydata) #yerr / ydata  #as long as it is constant this doesnt matter
        else:
            yerr = 1./weight
            
        #define fit function
        fitfunc = lambda p, x: p[0] - p[0] * np.exp(-p[1]*x)
        errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
        #guess fit results first
        pinit = [1.2, 1e3] #,1e3]

    else: 
        print "function: " + function + " not implemented"
        sys.exit(1)
        
    if method=='leastsq':
        #from IPython.core.debugger import Tracer ; Tracer()()
        out = optimize.leastsq(errfunc, pinit,
                        args=(x, y, yerr), full_output=1)
    elif method=='leastsq2D':
        #print x1.shape,x2.shape,y.shape,pinit
        #print x1,x2,y,pinit
        out = optimize.leastsq(errfunc, pinit,
                        args=(x, x2, y, yerr), full_output=1)
        #from IPython.core.debugger import Tracer ; Tracer()()
        
        #out = optimize.leastsq(errfunc, pinit,args=(x, y, yerr), full_output=1)
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
    elif func=='mod_powerlaw': #convert the fitted variables back 
        #pfinal[ = pfinal[1]
        pfinal[0] = 10.0**pfinal[0]

    elif func=='Atlas':
        pfinal = np.insert(pfinal,1,pfinal[0])

    
    return pfinal,covar # amp,index


def fitting_wrapper(ax,particle_dic,prop,N_mono_list,usecolors,fit_dic,diam,function="powerlaw",linewidthinvalid=0.05,linewidthvalid=2.0,allagg=True,noplotting=False):
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
        - diam: the diameter array on which the fit is analyzed
        - function: the functional relationship to which the property is fitted (e.g. powerlaw (default), Atlastype)
        - linewidthinvalid: linewidth at diameter ranges where there is only sparse data (particle from the aggregate model)
        - linewidthvalid: linewidth at diameter ranges where there are sufficient particle from the aggregate model
        - allagg: fit a line also for all the selected aggregates
        - noplotting: skip the plotting and calculate only the coefficients
    OUTPUT:
        - ax: the axis on which the fits are added as lines
    
    '''
    ####
    #fit a power-law for each quantitity and Atlas-type function for "vterm_" 
    ####


    #define some overall properties of the fit-dictionary and initialize it
    #set up a diameter array (for fitting indirectly (over the m-D,A-D coefficients) and plotting)
    #diam = np.logspace(-4,-1,100)
    fit_dic["diam"] = diam
    fit_dic["N_mono_list"] = N_mono_list

    for i_N_mono,N_mono_now in enumerate(N_mono_list): #loop over all monomer numbers
        if not any(particle_dic["N_monomer"]==N_mono_now): #if there are no particles with this monomer number just skip this loop
            break
        
        # make a fit for the current number of monomers (power-law; always directly to the points coming from the aggregate model)
        [fitresult,covar] = fit_data(particle_dic["diam"][particle_dic["N_monomer"]==N_mono_now],particle_dic[prop][particle_dic["N_monomer"]==N_mono_now],func=function,weight=fit_dic["dens" + "_Nmono_" + str(N_mono_now)])
            
        #save fit coeffients and calculate arrays of fitted masses and areas
        fit_dic[prop + "_coeff_Nmono_" + str(N_mono_now)] = fitresult #copy the A-D and m-D fit coefficients in the fit_dic dictionary
        fit_dic[prop + "_Nmono_" + str(N_mono_now)] = fit_dic[prop + "_coeff_Nmono_" + str(N_mono_now)][0]*fit_dic["diam"]**fit_dic[prop + "_coeff_Nmono_" + str(N_mono_now)][1] #calculate arrays of masses and areas 

        ###
        #plot the power-law (of mass or area)
        ###
        #seperate the plot according to the separation of sparse and dense enough data
        if not noplotting:
            ax.plot(fit_dic["diam"][fit_dic["diam_dens_enough_" + str(N_mono_now)]],fit_dic[prop + "_Nmono_" + str(N_mono_now)][fit_dic["diam_dens_enough_" + str(N_mono_now)]],color=usecolors[i_N_mono],linewidth=linewidthvalid) #here only the "valid" diameter ranges are plotted
            ax.plot(fit_dic["diam"],fit_dic[prop + "_Nmono_" + str(N_mono_now)],color=usecolors[i_N_mono],linewidth=linewidthinvalid)#the thin line can be plotted on the whole range

    if allagg==True:
        ##
        # make a fit for all aggregates in N_mono_list except N_mono=1 (power-law; always directly to the points coming from the aggregate model)
        ##
        N_mono_in_list = [False]*particle_dic["N_monomer"].shape[0]
        for N_mono_now in N_mono_list[1:]:
            N_mono_matching_now = (N_mono_now==particle_dic["N_monomer"])
            N_mono_in_list = np.logical_or(N_mono_matching_now,N_mono_in_list)
        #fit to m-D and A-D coefficients for all particles with N_monomer>1
        [fitresult,covar] = fit_data(particle_dic["diam"][N_mono_in_list],particle_dic[prop][N_mono_in_list],func='powerlaw')
        #save fit coeffients and calculate arrays of fitted masses and areas (for all particles with N_monomer>1)
        fit_dic[prop + "_coeff_Nmono_allagg"] = fitresult
        fit_dic[prop + "_Nmono_allagg"] = fit_dic[prop + "_coeff_Nmono_allagg"][0]*fit_dic["diam"]**fit_dic[prop + "_coeff_Nmono_allagg"][1] #calculate arrays of masses and areas based on the m-D fits
        if not noplotting:
            #plot the power-law (for all particles with N_monomer>1)
            ax.plot(fit_dic["diam"][fit_dic["diam_dens_enough_allagg"]],fit_dic[prop + "_Nmono_allagg"][fit_dic["diam_dens_enough_allagg"]],color='black',linewidth=linewidthvalid)
            fitline_allagg, = ax.plot(fit_dic["diam"],fit_dic[prop + "_Nmono_allagg"],color='black',linewidth=linewidthinvalid,label="Nmono>1")
        else:
            fitline_allagg=0 #dummy
    else:
        fitline_allagg=0 #dummy
        
    return ax, particle_dic,fitline_allagg #particle_dic is returned, because vterm_..._calculated via mD_AD is added

    
    
def add_fitresult_text(ax,fit_dic,prop,N_mono_list,function="powerlaw",hide_Nmono_fits=False,allagg=True,N_mono_listforlabel=None):
    '''
    adds the coefficients of the fit result in the form of a text
    INPUT:
        ax: the axis on which the text is added (this is also returned)
        fit_dic: a dictionary which contains the coefficients of the fits
        prop: name of the property (in the fit_dic dictionary) which is displayed now
        N_mono_list: list of the plotted monomer numbers
        function: function which was fitted
        hide_Nmono_fits: hide Nmono fits except Nmono=1 and Nmono>1 (can be useful when plotting many different monomer numbers)
        allagg: fit a line also for all the selected aggregates

    '''
    if N_mono_listforlabel==None: #use the standard N_mono_list for labelling (this is the case for non-grouped fits)
        N_mono_listforlabel=N_mono_list
    
    #initialize a string with the result of the powerlaw and  Atlas-type fits
    if function=="powerlaw":
        fitresult_string = 'fit a*D**b         a                 b\n' #string which contains the fit results
    elif function=="Atlas":
        fitresult_string =    'fit A-B*exp(-C*D) A                B            C\n' #string which contains the fit results
    for i_list,N_mono_now in enumerate(N_mono_list):
        if (not hide_Nmono_fits) or N_mono_now==1:
            fitresult_string += "Nmono=" + str(N_mono_listforlabel[i_list]) + ': ' + str(fit_dic[prop  + "_coeff_Nmono_" + str(N_mono_now)]) + '\n' #write the current fit result into the general string
    if allagg:
        fitresult_string += "Nmono>1: " + str(fit_dic[prop  + "_coeff_Nmono_allagg"]) + '\n'  #add also the result for all monomer numbers
    
    #add the text of the fitresult to the plot
    ax.text(0.0,1.0,fitresult_string,horizontalalignment='left',verticalalignment='top',transform = ax.transAxes,fontsize=6) #transform: give number in relative axes units instead of data coords

    
    return ax