'''
some outsourced functions for the processing of Jussis-aggregates
'''
import numpy as np
from numpy import pi, r_
import matplotlib.pyplot as plt
from scipy import optimize
import sys
from IPython.core.debugger import set_trace

def init_particle_dict():
    '''
    initialize the particle dictionary
    '''
    particle_dic = dict() #dictionary which contains information about the type of the pristine crystals (dendrites, needles, plates, ...) the underlying distribution of the pristine crystals and the properties of the aggregates (m,A,D,...)
    particle_dic["particle_type"] = [] #type of the pristine crystals
    particle_dic["N_monomer"] = [] #number of monomers
    particle_dic["mass"] = [] #particle mass
    particle_dic["area"] = [] #particle area
    particle_dic["diam"] = [] #particle diameter
    #particle_dic["sizeparam_index"] = [] #size parameter
    particle_dic["as_ratio"] = [] #aspect ratio
    particle_dic["area_partal10std"] = [] #particle area with a partial alignment (std-dev: 10deg)
    particle_dic["area_partal20std"] = [] #particle area with a partial alignment (std-dev: 20deg)
    particle_dic["area_partal40std"] = [] #particle area with a partial alignment (std-dev: 40deg)
    particle_dic["area_partal60std"] = [] #particle area with a partial alignment (std-dev: 60deg)
    particle_dic["N_initial"]        = [] #number of particles in initial distribution
    particle_dic["area_aligned_side_yz"]     = [] #side projected area of the aligned particle
    particle_dic["area_aligned_side_xz"]     = [] #side projected area of the aligned particle
    particle_dic["area_aligned_side_random"]     = [] #side projected area of the aligned particle
    particle_dic["area_aligned_filled"]     = [] #projected area of the aligned particle with filled holes
    particle_dic["Dmax_random_aligned"] = [] #maximum dimension from the "area_aligned_side_random"-random projection
    particle_dic["area_partal40_side"]     = [] #side projected area of the tumbled particle

    #particle_dic["HWJussi"] = [] #particle area
    #particle_dic["KCJussi"] = [] #particle diameter
    particle_dic["Dmax_xy"] = [] #maximum dimension projected to ..-.. plane
    particle_dic["Dmax_xz"] = [] #maximum dimension projected to ..-.. plane
    particle_dic["Dmax_yz"] = [] #maximum dimension projected to ..-.. plane
    particle_dic["Dmax_xy_aligned"] = [] #maximum dimension projected to ..-.. plane
    particle_dic["Dmax_xz_aligned"] = [] #maximum dimension projected to ..-.. plane
    particle_dic["Dmax_yz_aligned"] = [] #maximum dimension projected to ..-.. plane  
    particle_dic["rad_gyr"] = [] #radius of gyration
    particle_dic["D_circum_proj_area"] = [] #radius encircling the projected area
    return particle_dic

def calc_mD_AD_coeffs(mono_type,skip_mono=False):
    '''
    calculate the coefficients in the m-D relation (m=aD**b)
    INPUT: mono_type monomer habit name
    skip_mono: dont take the monomer and therefore ignore warnings
    OUTPUT: a_theor,b_theor: coefficients in m=a_theor*D**b_theor relation
    '''
    rho_i = 917.6 #define ice density
    '''
    if mono_type in ("dendrite","mixdendneedle"):
        #theoretical m-D relationships
        a_theor=0.096 #8.89e-2 #ATTENTION: this comes from the 1mum simulation fit
        b_theor=2.38 #ATTENTION: this comes from the 1mum simulation fit
        #theoretical A-D relationships
        c_theor=0.222 #ATTENTION: this comes from the 1mum simulation fit
        d_theor=2. #ATTENTION: this comes from the 1mum simulation fit
    elif mono_type=="plate": #copied here from riming.py
        #theoretical m-D relationships
        a_theor=rho_i*3.*np.sqrt(3)/2./4*(1.737e-3)*(1./2.)**0.474
        b_theor=2.+0.474
        #theoretical A-D relationships
        c_theor=0.675 #3.*np.sqrt(3)/2. *1./4. #1/4 comes from D=a/2 with a as the side of the hexagon
        d_theor=2.
    elif mono_type in ("needle"): #,"mixofdentneed"):
        #theoretical m-D relationships
        a_theor=rho_i*3.*np.sqrt(3)/2.*(1.319e-3)**2
        b_theor=2.*0.437+1
        #theoretical A-D relationships
        c_theor=0.0029#(1.+2./np.sqrt(2.))*1.319e-3 #2.4e-3 #
        d_theor=1.47 #1+0.437 #1.44
    elif mono_type=="rosette":
        #theoretical m-D relationships
        a_theor=0.7762081
        b_theor=2.48828817
        #theoretical A-D relationships
        c_theor=0.0412075
        d_theor=1.73699001
    elif mono_type=="bullet":
        
        #From Hong, 2007
        #alpha = 28*(np.pi/180.0)
        #f = np.sqrt(3)*1.552/np.tan(alpha)        
        #def L_func(L):
        #    return 2*L + f*L**0.63 - D*1e6
        #solve L from D numerically
        #self.L = brentq(L_func, 0, D*1e6/2.0) * 1e-6        
        #self.a = (1.552 * (self.L*1e6)**0.63) * 1e-6
        #self.t = np.sqrt(3)*self.a/(2*np.tan(alpha))
        #self.D = D
        
        #theoretical m-D relationships
        a_theor=np.nan
        b_theor=np.nan
        #theoretical A-D relationships
        c_theor=np.nan
        d_theor=np.nan
        if not skip_mono:
            raw_input("attention: the monomer power-laws are taken from columns!")
        #theoretical m-D relationships
        a_theor=0.07819 #rho_i*3.*np.sqrt(3)/2.*3.48**2*1e-6
        b_theor=2.14 #2.
        #theoretical A-D relationships
        c_theor=0.0143 #9.4e-3 #(1.+2./np.sqrt(2.))*3.48*1e-6**0.5
        d_theor=1.614 ##1.54 #1.5
    elif mono_type in ("column","mixcolumndend","mixcoldend1"):
        #theoretical m-D relationships
        a_theor=rho_i*3.*np.sqrt(3)/2.*3.48**2*1e-6
        b_theor=2.
        #theoretical A-D relationships
        c_theor=0.010 #(1.+2./np.sqrt(2.))*3.48*1e-6**0.5
        d_theor=1.58 #1.5 #a*L=L^0.5*L
    elif mono_type=="spheroid":
        #theoretical m-D relationships
        a_theor=np.pi/6*rho_i
        b_theor=3
        #theoretical A-D relationships
        c_theor=np.pi/4
        d_theor=2
    else:
        print mono_type + 'not defined in tools_for_processing_Jagg.calc_mD_AD_coeffs'
        sys.exit()
    print mono_type,a_theor,b_theor,c_theor,d_theor
    return a_theor,b_theor,c_theor,d_theor
    '''
    #take all from the 5mum run
    if mono_type in ("dendrite","mixdendneedle"):
        #m-D relationships
        a_theor=0.07368319
        b_theor=2.33289694
        #A-D relationships
        c_theor=0.14248231
        d_theor=1.93681387
    elif mono_type=="plate": #copied here from riming.py
        #m-D relationships
        a_theor=0.78784625
        b_theor=2.48057059
        #A-D relationships
        c_theor=0.63132449
        d_theor=1.99485988
    elif mono_type in ("needle"): #,"mixofdentneed"):
        #m-D relationships
        a_theor=0.00467768
        b_theor=1.89199142
        #A-D relationships
        c_theor=0.00220797
        d_theor=1.42440424
    elif mono_type=="rosette":
        #m-D relationships
        a_theor=0.54662064
        b_theor=2.43787266
        #A-D relationships
        c_theor=0.04434447
        d_theor=1.7480942
    elif mono_type=="bullet":
        print mono_type + 'not defined in tools_for_processing_Jagg.calc_mD_AD_coeffs'
        sys.exit()
    elif mono_type in ("column","mixcolumndend","mixcoldend1"):
        #m-D relationships
        a_theor=0.0458012
        b_theor=2.07368716
        #A-D relationships
        c_theor=0.00798099
        d_theor=1.53934525
    elif mono_type=="spheroid":
        a_theor=np.pi/6.*rho_i
        b_theor=3.
        c_theor=np.pi/2.
        d_theor=2.
    else:
        a_theor = 1e-100
        b_theor= 1e-100
        c_theor= 1e-100
        d_theor= 1e-100
        print mono_type + 'not defined in tools_for_processing_Jagg.calc_mD_AD_coeffs'
        #sys.exit()
    print mono_type,a_theor,b_theor,c_theor,d_theor
    return a_theor,b_theor,c_theor,d_theor

def read_particle_prop_files_rimed(prop_file_folder,particle_type):
    '''
    read the particle property files and save them in the "particle_dic" dictionary
    INPUT:  
            prop_file_folder: absolute path of the aggregates property file
            particle_type: habit
    OUTPUT: particle_dic: dictionary containing the properties of the individual simulated particles
    '''
    import os
    import glob
    import re
    import csv
    import sys
    particle_dic = init_particle_dict() #initialize dictionary which contains information about the type of the pristine crystals (dentrites, needles, plates, ...) the underlying distribution of the   
    print "processing " + particle_type 

    if not os.path.exists(prop_file_folder): #check if folder exists
        print prop_file_folder, " doesnt exist (exit in __tools_for_processing_Jagg.py"
        sys.exit(1)
        njk
    if len(glob.glob(prop_file_folder +  particle_type + '*properties.txt'))==0: #check if any file of this habit is there
        print prop_file_folder + "*", " doesnt exist (exit in __tools_for_processing_Jagg.py"
        sys.exit(1)
        fby 
    i_file=0
    for filename in glob.glob(prop_file_folder +  particle_type + '*properties.txt'): #loop over all property files (with different mean sizes)
        #verbose reading of files
        print "reading: " + filename
        with open(filename,"rb") as txtfile: #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in
            prop_reader = csv.reader(filter(lambda row: row[0]!='D',txtfile), delimiter=' ', quoting=csv.QUOTE_NONNUMERIC, lineterminator=os.linesep) #row[0]!='#': skip the header; quoting avoids '' for formatted string; lineterminator avoids problems with system dependend lineending format https://unix.stackexchange.com/questions/309154/strings-are-missing-after-concatenating-two-or-more-variable-string-in-bash?answertab=active#tab-top
            particle_dic["particle_type"] = np.append(particle_dic["particle_type"],particle_type)
            for i_row,row_content in enumerate(prop_reader):
                particle_dic["mass"] = np.append(particle_dic["mass"],row_content[1])

                particle_dic["diam"] = np.append(particle_dic["diam"],row_content[0])

                particle_dic["area"] = np.append(particle_dic["area"],row_content[2]) 
                particle_dic["N_monomer"] = np.append(particle_dic["N_monomer"],2.0) #dummy: all are aggregates 
    if False: 
        print "no particles found"
        sys.exit(1)                    
    N_mono_list = np.append(np.append(np.array(range(1,10)),np.array(range(10,100,10))),np.array(range(100,1001,100))) #the full monomer number list
    return particle_dic,N_mono_list



def read_particle_prop_files(prop_file_folder=  '/data/optimice/Jussis_aggregates_bugfixedrotation/',
                                                D_small_notuse=1e-4, #ATTENTION: particle smaller than D_small_notuse are not considered (e.g. because of resolution issues)
                                                N_small_notuse=1,    #ATTENTION:  monomer numbers smaller than N_small_notuse are not considered
                                                grid_res_array = [1e-6,5e-6,10e-6], #array of grid resolutions (defines the files which are read in)
                                                particle_type = "plate", #define the habit
                                                test_with_small_sample = False, #read only 1/10 of the data in for faster plotting
                                                use_only_monomers = False, #get only Nmono=1 particle properties back
                                                forKamil = False
                                                ):
    '''
    read the particle property files and save them in the "particle_dic" dictionary
    INPUT:  
            prop_file_folder: absolute path of the aggregates property file
            D_small_notuse: particles smaller than that value (in m) are not used
            N_small_notuse: monomer numbers smaller than N_small_notuse are not considered
            grid_res_array: array with grid resolution
            use_only_monomers :get only Nmono=1 particle properties back
    
            forKamil: some different format for Kamils closure paper
    OUTPUT: particle_dic: dictionary containing the properties of the individual simulated particles
    '''
    import os
    import glob
    import re
    import csv
    import sys
    particle_dic = init_particle_dict() #initialize dictionary which contains information about the type of the pristine crystals (dentrites, needles, plates, ...) the underlying distribution of the   
    for grid_res in grid_res_array: #loop over different grid-resolutions in order to have small monomer numbers modelled woth higher resolution
        print "processing " + particle_type + " with resolution: ",grid_res
        sensrun_folder = 'res_'+ str(int(grid_res*1e6)) + 'mum/'
        if not os.path.exists(prop_file_folder + sensrun_folder): #check if folder exists
            print prop_file_folder + sensrun_folder, " doesnt exist (exit in __tools_for_processing_Jagg.py"
            sys.exit()
        if len(glob.glob(prop_file_folder + sensrun_folder +  particle_type + '*properties.txt'))==0: #check if any file of this habit is there
            print prop_file_folder + sensrun_folder + particle_type + "*", " doesnt exist (exit in process_Jussis_aggregate_model_mod_powerlaw.py"
            sys.exit()
        
        if test_with_small_sample:
            test_subsample = "*1E"
        else:
            test_subsample = ""
        i_file=0
        for filename in glob.glob(prop_file_folder + sensrun_folder +  particle_type + test_subsample + '*properties.txt'): #loop over all property files (with different mean sizes)
            i_file+=1
            #read size parameter from filename in order to be able to disregard some size parameter
            m=re.search(prop_file_folder + sensrun_folder + particle_type + '_(.+?)_properties.txt', filename) #search for patter
            size_param_now = float(m.group(1))
            #if size_param_now<-999: #ignore all sizeparameter above ... or below ... range:[50,50..,51.,..,1000]
            #print size_param_now
            #if size_param_now<400 or size_param_now>600: #take only very narrow sizeparameter
            #    continue
        
            #verbose reading of files
            print "reading: " + filename
            with open(filename,"rb") as txtfile: #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in
                prop_reader = csv.reader(filter(lambda row: row[0]!='#',txtfile), delimiter=' ', quoting=csv.QUOTE_NONNUMERIC, lineterminator=os.linesep) #row[0]!='#': skip the header; quoting avoids '' for formatted string; lineterminator avoids problems with system dependend lineending format https://unix.stackexchange.com/questions/309154/strings-are-missing-after-concatenating-two-or-more-variable-string-in-bash?answertab=active#tab-top

                #define the monomer numbers to read in
                take_all_Nmono_bool = True
                if grid_res==1e-6: #particles < 10mum are modelled only for small monomer numbers
                    N_mono_list_tmp = np.array(range(2,10))
                elif grid_res==5e-6: #particles < 10mum are modelled only for small monomer numbers
                    if 1e-6 in grid_res_array:
                        N_mono_list_tmp = np.array(range(10,101,10))
                    else:
                        N_mono_list_tmp = np.array(np.append(range(1,11,1),np.array(range(10,101,10))))
                elif grid_res==10e-6:#take only the large monomer numbers from the low-resolution runs
                    if 5e-6 in grid_res_array:
                        N_mono_list_tmp = np.array(range(100,1001,100))
                    else:
                        if 1e-6 in grid_res_array:
                            N_mono_list_tmp = np.array(np.append(np.array(range(10,101,10)),np.array(range(100,1001,100))))
                        else:
                            N_mono_list_tmp = np.array(np.append(np.array(range(1,11,1)),np.append(np.array(range(10,101,10)),np.array(range(100,1001,100)))))
                            #N_mono_list_tmp = np.array([1])
                else: #some tests with higher/lower resolution
                    N_mono_list_tmp = np.array(np.append(np.append(range(1,11,1),np.array(range(10,101,10))),range(200,1001,100)))
                for i_row,row_content in enumerate(prop_reader): #TODO: the row is not any more equal to the monomer number #read the row with N_mono_now (the currently considered monomer number)
                    if forKamil and len(row_content)<6:
                        pass

                    ################
                    ##limit N_mono##
                    ################
                    #if row_content[0]>20 and (row_content[0] in N_mono_list_tmp): #i_row==0 or i_row==4 or i_row==9: #limit N_mono 
                    if (row_content[0] in N_mono_list_tmp): #i_row==0 or i_row==4 or i_row==9:
                        if (row_content[3]<D_small_notuse): #check size limit
                            continue
                        if (row_content[0]<N_small_notuse): #check Nmono limit
                            continue
                        if use_only_monomers:
                            if (row_content[0]>1): #only monomers
                                continue
                            if (row_content[3]>1e-3): #only small particles 
                                continue

                        particle_dic["particle_type"] = np.append(particle_dic["particle_type"],particle_type)
                        particle_dic["N_monomer"] = np.append(particle_dic["N_monomer"],row_content[0]) #python indices start with 0
                        particle_dic["mass"] = np.append(particle_dic["mass"],row_content[1])

                        particle_dic["diam"] = np.append(particle_dic["diam"],row_content[3]) #np.sqrt(row_content[4]/np.pi)) #ATTENTION: test take area equivalent diameter

                        if forKamil:
                            particle_dic["area"] = np.append(particle_dic["area"],row_content[2]) #,row_content[2]) -[2] for tumbling [4] for disabled tumbling
                        else:

                            particle_dic["area"] = np.append(particle_dic["area"],row_content[4]) #,row_content[2]) -[2] for tumbling [4] for disabled tumbling
                        try:
                            particle_dic["as_ratio"] = np.append(particle_dic["as_ratio"],row_content[5]) 
                        except:
                            pass
                        try:
                            particle_dic["area_partal10std"] = np.append(particle_dic["area_partal10std"],row_content[6]) 
                        except:
                            pass
                        try:
                            particle_dic["area_partal20std"] = np.append(particle_dic["area_partal20std"],row_content[7])
                        except:
                            pass
                        particle_dic["area_partal40std"] = np.append(particle_dic["area_partal40std"],row_content[2])
                        try:
                            particle_dic["area_partal60std"] = np.append(particle_dic["area_partal60std"],row_content[8])
                        except:
                            pass
                        try:
                            particle_dic["area_aligned_side_yz"] = np.append(particle_dic["area_aligned_side_yz"],row_content[9])
                        except:
                            pass
                        try:
                            particle_dic["area_partal40_side"] = np.append(particle_dic["area_partal40_side"],row_content[10])
                        except:
                            pass
                        #try to read some lines, that might not always be there
                        #try:
                        #    particle_dic["area_side_aligned"] = np.append(particle_dic["N_initial"],row_content[10])
                        #except:
                        #    pass
                        #try:
                        #    particle_dic["N_initial"] = np.append(particle_dic["N_initial"],row_content[11])
                        #except:
                        #    pass
                        #try:
                        #    particle_dic["N_initial"] = np.append(particle_dic["N_initial"],row_content[11])
                        #except:
                        #    pass
                        try:
                            particle_dic["Dmax_xy"] = np.append(particle_dic["Dmax_xy"],row_content[11])
                        except:
                            pass
                        try:
                            particle_dic["Dmax_xz"] = np.append(particle_dic["Dmax_xz"],row_content[12])
                        except:
                            pass
                        try:
                            particle_dic["Dmax_yz"] = np.append(particle_dic["Dmax_yz"],row_content[13])
                        except:
                            pass
                        try:
                            particle_dic["Dmax_xy_aligned"] = np.append(particle_dic["Dmax_xy_aligned"],row_content[14])
                        except:
                            pass
                        try:
                            particle_dic["Dmax_xz_aligned"] = np.append(particle_dic["Dmax_xz_aligned"],row_content[15])
                        except:
                            pass
                        try:
                            particle_dic["Dmax_yz_aligned"] = np.append(particle_dic["Dmax_yz_aligned"],row_content[16])
                        except:
                            pass
                        try:
                            particle_dic["rad_gyr"] = np.append(particle_dic["rad_gyr"],row_content[17])
                        except:
                            pass

                        try:
                            particle_dic["D_circum_proj_area"] = np.append(particle_dic["D_circum_proj_area"],row_content[18])
                        except:
                            pass
                        try:
                            particle_dic["area_aligned_side_xz"] = np.append(particle_dic["area_aligned_side_xz"],row_content[19])
                        except:
                            pass
                        try:
                            particle_dic["area_aligned_side_random"] = np.append(particle_dic["area_aligned_side_random"],row_content[20])
                        except:
                            pass
                        try:
                            particle_dic["Dmax_random_aligned"] = np.append(particle_dic["Dmax_random_aligned"],row_content[21])
                        except:
                            pass
                        try:
                            particle_dic["area_aligned_filled"] = np.append(particle_dic["area_aligned_filled"],row_content[22])
                        except:
                            pass
                N_mono_list = np.append(np.append(np.array(range(1,10)),np.array(range(10,100,10))),np.array(range(100,1001,100))) #the full monomer number list
        if False: #len(particle_dic["diam"])==0:
            print "no particles found"
            sys.exit()                    
    return particle_dic,N_mono_list


def fit_data(xdata,ydata,func='powerlaw',method='leastsq',weight='None',habit='None',prop='None'):
    '''source: https://scipy-cookbook.readthedocs.io/items/FittingData.html
    fit any 2-dimensional data with different methods and target functions
    xdata: ordinate points
    ydata: coordinate points
    func: select the function to which the data should be fitted
    habit: just used for modified power-law (get parameter from monomer)
    prop: just used for modified power-law (parameters are different for mass and area)
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
        if weight=="None":
            weight = 1./np.ones_like(ydata)
        else:
            pass #from IPython.core.debugger import Tracer ; Tracer()()
        yerr=1./(weight)
        #yerr=np.log10(yerr)
        #from IPython.core.debugger import Tracer ; Tracer()()

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
        if habit=='None':
            print "a habit should be defined for the mod_powerlaw fit"
            sys.exit()
        if prop=='None':
            print "the fitted property should be defined"
            sys.exit()


        #get m-D from the assumptions in the aggregate model
        a,b,c,d = calc_mD_AD_coeffs(habit)
        
        #from IPython.core.debugger import Tracer ; Tracer()()
        if prop=="mass":
            x = np.log10(xdata[0]) #diameter
            x2 = (xdata[1]) #monomer number
            #first normalize y by the properties of the monomer
            ydata = ydata/(a*xdata[0]**b)
            #from IPython.core.debugger import Tracer ; Tracer()()
            # define our (line) fitting function #inspired by: https://stackoverflow.com/questions/29374296/two-dimensional-fit-with-python
            fitfunc = lambda p, x, x2: 1+p[0]*np.log10(x2) #(p[0]-p[0]*np.exp(p[1]*(x2-1)))+ p[2]*x  #+ p[1]*x2/(1+p[0]*(x2)) #*x**(p[1]*np.log10(x2)) * np.log10(a)*x**b #*(p[2])#p[0]*(1+p[1]*(x2-1)) + p[2]*(1+p[3]/(x2))* x
            errfunc = lambda p, x, x2, y, err: (y - fitfunc(p, x, x2)) / err
        elif prop=="area":
            x = np.log10(xdata[0]) #diameter
            x2 = np.log10(xdata[1]) #monomer number
            #first normalize y by the properties of the monomer
            ydata = ydata/(a*xdata[0]**b)
            # define our (line) fitting function #inspired by: https://stackoverflow.com/questions/29374296/two-dimensional-fit-with-python
            fitfunc = lambda p, x, x2: (1+p[0]*(x2+1)) + np.log10(c) + d*x # + p[1]*x2*x #*x**(p[1]*np.log10(x2)) * np.log10(a)*x**b #*(p[2])#p[0]*(1+p[1]*(x2-1)) + p[2]*(1+p[3]/(x2))* x
            errfunc = lambda p, x, x2, y, err: (y - fitfunc(p, x, x2)) / err

        if weight=="None":
            weight = np.ones_like(ydata)   #ATTENTION: weighting by McSnow histogram turned off

        #define weight
        yerr=1./(ydata*weight)
        yerr = np.log10(yerr)
        #logarithmate the y-data
        y = np.log10(ydata)
            
        #guess fit results first
        pinit = [2.0,-1e-3,0] #[np.log10(0.01), 0.2,2.5,-0.2]


    elif func=='Atlas':
        #Atlas-type fitting
        # y=A-B*exp(-gam*x)
        
        #stay linear for fitting
        x = (xdata)
        y = (ydata)
        
        #from IPython.core.debugger import Tracer ; Tracer()()
        if weight=='None': #dont apply any weighting
            yerr = np.ones_like(ydata) #yerr / ydata  #as long as it is constant this doesnt matter
        else:
            yerr = 1./weight
            
        #define fit function
        fitfunc = lambda p, x: p[0] - p[1] * np.exp(-p[2]*x)
        #fitfunc = lambda p, x: p[0] - p[0] * np.exp(-p[1]*x) #if two parameters only

        errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
        #guess fit results first
        pinit = [1.2, 1.2 ,1e3]
        #pinit = [1.2, 1e3] #if two parameters only


    else: 
        print "function: " + function + " not implemented"
        sys.exit(1)
        
    if method=='leastsq':
        #from IPython.core.debugger import Tracer ; Tracer()()
        out = optimize.leastsq(errfunc, pinit,
                        args=(x, y, yerr), full_output=1)
    elif method=='leastsq2D':

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
        pass #pfinal[0] = 10.0**pfinal[0]

    elif func=='Atlas':
        pass #
        #pfinal = np.insert(pfinal,1,pfinal[0]) #if two parameters only

    
    return pfinal,covar # amp,index


def fitting_wrapper(ax,particle_dic,prop,N_mono_list,usecolors,fit_dic,diam,function="powerlaw",linewidthinvalid=0.9,linewidthvalid=2.0,allagg=True,plot_fits=False,plot_binary_fits=False, weight='None'):
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
        - plot_fits: plot also the fit line for all monomer numbers and calculate the coefficients for all monomer numbers
        - plot_binary_fits: plot fits for the monomer and all aggregates but not necessarily for all monomer numbers (see plot_fits)
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
            #raw_input("there is none of " + str(N_mono_now))
            continue
        
        # make a fit for the current number of monomers (power-law; always directly to the points coming from the aggregate model)
        if not isinstance(weight,str) and weight=='None':
            [fitresult,covar] = fit_data(particle_dic["diam"][particle_dic["N_monomer"]==N_mono_now],particle_dic[prop][particle_dic["N_monomer"]==N_mono_now],func=function)
            print "check if you really want to get here 1 ( in tools_for_processing_Jagg.fitting_wrapper";from IPython.core.debugger import Tracer ; Tracer()()
        else :#apply a weight
            [fitresult,covar] = fit_data(particle_dic["diam"][particle_dic["N_monomer"]==N_mono_now],particle_dic[prop][particle_dic["N_monomer"]==N_mono_now],func=function,weight=weight)
            
        '''
        else: 
            [fitresult,covar] = fit_data(particle_dic["diam"][particle_dic["N_monomer"]==N_mono_now],particle_dic[prop][particle_dic["N_monomer"]==N_mono_now],func=function)
            print "check if you really want to get here 2 ( in tools_for_processing_Jagg.fitting_wrapper";from IPython.core.debugger import Tracer ; Tracer()()
        '''
            
        #save fit coeffients and calculate arrays of fitted masses and areas
        fit_dic[prop + "_coeff_Nmono_" + str(N_mono_now)] = fitresult #copy the A-D and m-D fit coefficients in the fit_dic dictionary
        fit_dic[prop + "_Nmono_" + str(N_mono_now)] = fit_dic[prop + "_coeff_Nmono_" + str(N_mono_now)][0]*fit_dic["diam"]**fit_dic[prop + "_coeff_Nmono_" + str(N_mono_now)][1] #calculate arrays of masses and areas 

        ###
        #plot the power-law (of mass or area)
        ###
        #seperate the plot according to the separation of sparse and dense enough data
        if plot_binary_fits and N_mono_now==1:
            #ax.plot(fit_dic["diam"][fit_dic["diam_dens_enough_" + str(N_mono_now)]],fit_dic[prop + "_Nmono_" + str(N_mono_now)][fit_dic["diam_dens_enough_" + str(N_mono_now)]],color=usecolors[i_N_mono],linewidth=linewidthvalid) #here only the "valid" diameter ranges are plotted
            ax.plot(fit_dic["diam"],fit_dic[prop + "_Nmono_" + str(N_mono_now)],color=usecolors[i_N_mono],linewidth=linewidthinvalid,label="Nmono=1") #the thin line can be plotted on the whole range
            
    if allagg==True:
        ##
        # make a fit for all aggregates in N_mono_list except N_mono=1 (power-law; always directly to the points coming from the aggregate model)
        ##
        N_mono_in_list = [False]*particle_dic["N_monomer"].shape[0]
        for N_mono_now in N_mono_list[1:]:
            N_mono_matching_now = (N_mono_now==particle_dic["N_monomer"])
            N_mono_in_list = np.logical_or(N_mono_matching_now,N_mono_in_list)
        #fit to m-D and A-D coefficients for all particles with N_monomer>1
        [fitresult,covar] = fit_data(particle_dic["diam"][N_mono_in_list],particle_dic[prop][N_mono_in_list],func='powerlaw',weight=weight) 
        fit_dic[prop + "_coeff_Nmono_allagg"] = fitresult
        fit_dic[prop + "_Nmono_allagg"] = fit_dic[prop + "_coeff_Nmono_allagg"][0]*fit_dic["diam"]**fit_dic[prop + "_coeff_Nmono_allagg"][1] #calculate arrays of masses and areas based on the m-D fits
        if plot_binary_fits:
            #plot the power-law (for all particles with N_monomer>1)
            #ax.plot(fit_dic["diam"][fit_dic["diam_dens_enough_allagg"]],fit_dic[prop + "_Nmono_allagg"][fit_dic["diam_dens_enough_allagg"]],color='black',linewidth=linewidthvalid)
            fitline_allagg, = ax.plot(fit_dic["diam"],fit_dic[prop + "_Nmono_allagg"],color='y',linewidth=linewidthinvalid,label="Nmono>1")
        else:
            fitline_allagg=0 #dummy
    else:
        fitline_allagg=0 #dummy
        
    return ax, particle_dic,fitline_allagg #particle_dic is returned, because vterm_..._calculated via mD_AD is added



'''
START: scripts of the "polynomial/rational function powerlaw" fit
'''
def rational(x, p, q):
    """
    The general rational function description.
    p is a list with the polynomial coefficients in the numerator
    q is a list with the polynomial coefficients (except the first one)
    in the denominator
    The zeroth order coefficient of the denominator polynomial is fixed at 1.
    Numpy stores coefficients in [x**2 + x + 1] order, so the fixed
    zeroth order denominator coefficent must comes last. (Edited.)
    """
    return np.polynomial.polynomial.polyval(x, p) / np.polynomial.polynomial.polyval(x,q)


#polynom

def polynom1_powerlaw(x, p0, p1, r0, r1):
    '''
    simple powerlaw with a and b in (a*x**b) depending on a second parameter (e.g. the monomer number)
    '''
    x1,x2 = x
 
    return (p0+(p1*x2))+(r0+(r1*x2))*x1

def polynom1_powerlaw_fixp0_fixr0(x, p1, r1):
    '''
    simple powerlaw with a and b in (a*x**b) depending on a second parameter (e.g. the monomer number)
    '''
    x1,x2 = x

    return (p1*x2)+(r1*x2)*x1


def polynom2_powerlaw_fixp00_and_fixq00(x, p1, p2, r1, r2):

    x1,x2 = x

    #from IPython.core.debugger import Tracer ; Tracer()()
    return (1+p1*x2 + p2*x2**2)+(r1*x2 + r2*x2**2)*x1

#rational

def powerlaw_rational1_fixp0_fixr0(x, p1, q1, r1, s1):

    x1,x2 = x

    return rational(x2, [0.0, p1], [1.0, q1]) + rational(x2, [0.0, r1], [1.0, s1])*x1
    #return rational(x2, [0.0, 1.0], [1.0, 0.0]) + rational(x2, [0.0, 1.0], [1.0, 0.0])*x1

def powerlaw_rational1(x, p0, p1, q1, r0, r1, s1):
    
    x1,x2 = x
 
    return rational(x2, [p0, p1], [1.0, q1]) + rational(x2, [r0, r1], [1.0, s1])*x1

def powerlaw_rational2_fixp0_fixr0(x, p1, p2, q1, q2, r1, r2,  s1, s2):
    
    x1,x2 = x
 
    return rational(x2, [0.0, p1, p2], [1.0, q1, q2]) + rational(x2, [0.0, r1, r2], [1.0, s1, s2])*x1

def powerlaw_rational2(x, p0, p1, p2, q1, q2, r0, r1, r2,  s1, s2):
    
    x1,x2 = x
 
    return rational(x2, [p0, p1, p2], [1.0, q1, q2]) + rational(x2, [r0, r1, r2], [1.0, s1, s2])*x1

def powerlaw_rational3_fixp0_fixr0(x, p1, p2, p3, q1, q2, q3, r1, r2, r3,  s1, s2, s3):
    
    x1,x2 = x
 
    return rational(x2, [0.0, p1, p2, p3], [1.0, q1, q2, q3]) + rational(x2, [0.0, r1, r2, r3], [1.0, s1, s2, s3])*x1

def powerlaw_rational3(x,  p0, p1, p2, p3, q1, q2, q3, r0, r1, r2,r3,  s1, s2, s3):
    
    x1,x2 = x
 
    return rational(x2, [p0, p1, p2, p3], [1.0, q1, q2, q3]) + rational(x2, [r0, r1, r2, r3], [1.0, s1, s2, s3])*x1


def fit_2D_rational_powerlaw_function(xdata,ydata,func='rational1_1_powerlaw_fixp00_and_fixq00',method='lm',weight='None',habit='None',prop='None',fixp00=False,fixq00=False,guess=None,ftol=1e-8):

    ''' fit a rational function ( https://en.wikipedia.org/wiki/Polynomial_and_rational_function_modeling ) for a 2-dimensional dataset
        #code-parts are taken from:
            https://stackoverflow.com/questions/28372597/python-curve-fit-with-multiple-independent-variables
            https://stackoverflow.com/questions/29815094/rational-function-curve-fitting-in-python
        xdata: ordinate points (2-dimensional)
        ydata: coordinate points
        method='lm': Levenberg-Marquardt algorithm through leastsq (as in Frick. et. al 2010, GMD)
        func: select the function to which the data should be fitted (determining e.g. the order of the polynom)
        habit: used to normalize by the monomer habit (get parameter from monomer)
        prop: used to normalize by the monomer habit (parameters are different for mass and area)
    '''    
    
    if func=="powerlaw_polynom1_fixp0_fixr0":

        fitcoeffs, pcov = optimize.curve_fit(powerlaw_polynom1_fixp0_fixr0, xdata, ydata,method=method, p0=(0.0 , 0.0))   #[a, b]    
        
        print "\n\n\n\nfitcoeffs: \n", "     p1,           p2, \n", fitcoeffs
    
    elif func=="powerlaw_polynom1":

        fitcoeffs, pcov = optimize.curve_fit(powerlaw_polynom1, xdata, ydata,method=method, p0=(0.0 , 0.0, 0.0, 0.0))   #[a, b]    
        
        print "\n\n\n\nfitcoeffs: \n", "     p1,           p2, \n", fitcoeffs
    
    elif func=="powerlaw_rational1_fixp0_fixr0":

        if guess==None:
            fitcoeffs, pcov = optimize.curve_fit(powerlaw_rational1_fixp0_fixr0, xdata, ydata,method=method, p0=( 0.0,   #[p1]
                                                                                                            0.0,   #[q1]
                                                                                                            0.0,   #[r1]
                                                                                                            0.0,   #[s1]
                                                                                                            ))
        else:
            #print ftol,guess; raw_input()
            fitcoeffs, pcov = optimize.curve_fit(powerlaw_rational1_fixp0_fixr0, xdata, ydata,ftol=ftol,method=method, p0=( guess[0],   #[p1]
                                                                                                            guess[1],   #[q1]
                                                                                                            guess[2],   #[r1]
                                                                                                            guess[3],   #[s1]
                                                                                                            ))
        print "\n\n\n\nfitcoeffs: \n", "     p1,           q2,          r1,       s1  \n", fitcoeffs
    
    elif func=="powerlaw_rational1":

        fitcoeffs, pcov = optimize.curve_fit(powerlaw_rational1, xdata, ydata,method=method, p0=( 0.0,0.0,  #[p0,p1]
                                                                                                0.0,    #[q1]
                                                                                                0.0,0.0,    #[r0,r1]
                                                                                                0.0,    #[s0,s1]
                                                                                                        ))
        print "\n\n\n\nfitcoeffs: \n", fitcoeffs

    '''
    elif func=="powerlaw_rational2_fixp0_fixr0":

        #first perform the fit for one-degree polynomial less
        fitcoeffs_for_guess, pcov = optimize.curve_fit(powerlaw_rational1_fixp0_fixr0, xdata, ydata,method=method, p0=( 0.0,   #[p1]
                                                                                                        0.0,   #[q1]
                                                                                                        0.0,   #[r1]
                                                                                                        0.0,   #[s1]
                                                                                                        ))

        fitcoeffs, pcov = optimize.curve_fit(powerlaw_rational2_fixp0_fixr0, xdata, ydata,method=method, p0=( fitcoeffs_for_guess[0], 0.0,  #[p1,p2]
                                                                                                fitcoeffs_for_guess[1], 0.0,    #[q1,q2]
                                                                                                fitcoeffs_for_guess[2], 0.0,    #[r1,r2]
                                                                                                fitcoeffs_for_guess[3], 0.0,    #[s1,s2]
                                                                                                        ), maxfev=10000)
        
        print "\n\n\n\nfitcoeffs: \n", fitcoeffs

    elif func=="powerlaw_rational2":

        #first perform the fit for one-degree polynomial less
        fitcoeffs_for_guess, pcov = optimize.curve_fit(powerlaw_rational1, xdata, ydata,method=method, p0=( 0.0,0.0,  #[p0,p1]
                                                                                                0.0,    #[q1]
                                                                                                0.0,0.0,    #[r0,r1]
                                                                                                0.0,    #[s0,s1]
                                                                                                        ))

        fitcoeffs, pcov = optimize.curve_fit(powerlaw_rational2, xdata, ydata,method=method, p0=( 0.0,fitcoeffs_for_guess[0],fitcoeffs_for_guess[1],    #[p0,p1,p2]
                                                                                        fitcoeffs_for_guess[2], 0.0,   #[q1,q2]
                                                                                    0.0,fitcoeffs_for_guess[3], 0.0,   #[r0,r1,r2]
                                                                                        fitcoeffs_for_guess[4], 0.0,   #[s1,s2]
                                                                                                        ), maxfev=10000)
        print "\n\n\n\nfitcoeffs: \n", fitcoeffs

    elif func=="powerlaw_rational3_fixp0_fixr0":

        #first perform the fit for one-degree polynomial less
        fitcoeffs_for_guess, pcov = optimize.curve_fit(powerlaw_rational1_fixp0_fixr0, xdata, ydata,method=method, p0=( 0.0,   #[p1]
                                                                                                        0.0,   #[q1]
                                                                                                        0.0,   #[r1]
                                                                                                        0.0,   #[s1]
                                                                                                        ))
        
        fitcoeffs_for_guess2, pcov = optimize.curve_fit(powerlaw_rational2_fixp0_fixr0, xdata, ydata,method=method, p0=( fitcoeffs_for_guess[0], 0.0,    #[p0,p1,p2]
                                                                                fitcoeffs_for_guess[1], 0.0,   #[q1,q2]
                                                                                fitcoeffs_for_guess[2], 0.0,   #[r0,r1,r2]
                                                                                fitcoeffs_for_guess[3], 0.0,   #[s1,s2]
                                                                                                ), maxfev=10000)

        fitcoeffs, pcov = optimize.curve_fit(powerlaw_rational3_fixp0_fixr0, xdata, ydata,method=method, p0=( fitcoeffs_for_guess2[0], fitcoeffs_for_guess2[1], 0.0,    #[p1,p2,p3]
                                                                                                fitcoeffs_for_guess2[2], fitcoeffs_for_guess2[3], 0.0,    #[q1,q2,q3]
                                                                                                fitcoeffs_for_guess2[4], fitcoeffs_for_guess2[5], 0.0,    #[r1,r2,r3]
                                                                                                fitcoeffs_for_guess2[6], fitcoeffs_for_guess2[7], 0.0,    #[s1,s2,s3]
                                                                                                        ), maxfev=10000)
        
        print "\n\n\n\nfitcoeffs: \n", fitcoeffs

    elif func=="powerlaw_rational3":

        #first perform the fit for one-degree polynomial less
        fitcoeffs_for_guess, pcov = optimize.curve_fit(powerlaw_rational1, xdata, ydata,method=method, p0=( 0.0,0.0,  #[p0,p1]
                                                                                                0.0,    #[q1]
                                                                                                0.0,0.0,    #[r0,r1]
                                                                                                0.0,    #[s0,s1]
                                                                                                        ))
        print fitcoeffs_for_guess
        fitcoeffs_for_guess2, pcov = optimize.curve_fit(powerlaw_rational2, xdata, ydata,method=method, p0=( fitcoeffs_for_guess[0], fitcoeffs_for_guess[1], 0.0,  #[p0,p1,p2]
                                                                                                fitcoeffs_for_guess[2], 0.0,    #[q1,q2]
                                                                                                fitcoeffs_for_guess[3], fitcoeffs_for_guess[4], 0.0,    #[r0,r1,r2]
                                                                                                fitcoeffs_for_guess[5], 0.0,    #[s1,s2]
                                                                                                        ), maxfev=10000)
        
        print fitcoeffs_for_guess2

        fitcoeffs, pcov = optimize.curve_fit(powerlaw_rational3, xdata, ydata,method=method, p0=( fitcoeffs_for_guess2[0],fitcoeffs_for_guess2[1],fitcoeffs_for_guess2[2], 0.0,   #[p0,p1,p2,p3]
                                                                                                fitcoeffs_for_guess2[3],fitcoeffs_for_guess2[4], 0.0,    #[q1,q2,q3]
                                                                                                fitcoeffs_for_guess2[5],fitcoeffs_for_guess2[6],fitcoeffs_for_guess2[7], 0.0,    #[r0,r1,r2,r3]
                                                                                                fitcoeffs_for_guess2[8],fitcoeffs_for_guess2[9], 0.0,    #[s0,s1,s2,s3]
                                                                                                        ), maxfev=10000)
        print fitcoeffs

        print "\n\n\n\nfitcoeffs: \n", fitcoeffs
        '''
    
        
    return fitcoeffs


"""
END: scripts of the "rational function powerlaw" fit
"""

def fit_2D_rational_function_transformwrapper(diam,Nmono,ydata,func='powerlaw',method='lm',weight='None',habit='None',prop='None',guess='None',ftol=1e-8):
    """
    here we do the transformation of the variables and the "educated guess"
    """
    #get m-D from the assumptions in the aggregate model
    a,b,c,d = calc_mD_AD_coeffs(habit)
    
    #transform the data
    diam_log = np.log10(diam)
    Nmono_log = np.log10(Nmono)
    
    if prop=="mass":
        #first normalize y by the properties of the monomer
        ydata_normed = np.log10(ydata/(a*diam**b))
    elif prop=="area":
        #first normalize y by the properties of the monomer
        ydata_normed = np.log10(ydata/(c*diam**d))

    if "powerlaw" in func:
        fitcoeffs = fit_2D_rational_powerlaw_function(np.array([diam_log,Nmono_log]),ydata_normed,func=func,method='lm',weight=weight,habit=habit,prop=prop,guess=guess,ftol=ftol)
    #elif "rational" in func: #ATTENTION: THIS IS FOR FITTING A RATIONAL FUNCTION to m/m_mono (NOT CONSERVING THE POWER-LAW!), #functions for this method are currently commented
    #    fitcoeffs = fit_2D_rational_function(np.array([diam_log,Nmono_log]),ydata_normed,func=func,method='lm',weight=weight,habit=habit,prop=prop)

    return fitcoeffs

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
        fit_suffix = "" #dirty way to access right fit
    elif function=="powerlaw_fallspeed":
        fitresult_string = 'fit a*D**b         a                 b\n' #string which contains the fit results
        fit_suffix = "_powerlaw" #dirty way to access right fit
    elif function=="Atlas":
        fitresult_string =    'fit A-B*exp(-C*D) A                B            C\n' #string which contains the fit results
        fit_suffix = "" #dirty way to access right fit
    for i_list,N_mono_now in enumerate(N_mono_list):
        if (not hide_Nmono_fits) or N_mono_now==1:
            fitresult_string += "Nmono=" + str(N_mono_listforlabel[i_list]) + ': ' + str(fit_dic[prop  + "_coeff_Nmono_" + str(N_mono_now) + fit_suffix]) + '\n' #write the current fit result into the general string
    if allagg:
        fitresult_string += "Nmono>1: " + str(fit_dic[prop  + "_coeff_Nmono_allagg" + fit_suffix]) + '\n'  #add also the result for all monomer numbers
    
    if function=="powerlaw_fallspeed":
        #add the text of the fitresult to the plot
        ax.text(0.8,0.2,fitresult_string,horizontalalignment='left',verticalalignment='top',transform = ax.transAxes,fontsize=6) #transform: give number in relative axes units instead of data coords    
    if function=="Atlas" or function=="powerlaw":
        #add the text of the fitresult to the plot
        ax.text(0.0,1.0,fitresult_string,horizontalalignment='left',verticalalignment='top',transform = ax.transAxes,fontsize=6) #transform: give number in relative axes units instead of data coords

    
    return ax
