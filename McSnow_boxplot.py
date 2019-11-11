
import subprocess
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import __postprocess_McSnow
import __postprocess_SB
import McSnow_boxruns
from matplotlib import rc 
'''
this code executes McSnow_boxruns.py with different settings (particle properties,size distribution parameter, initialization (nrp0,iwc), supersaturation) and visualizes its output
'''
from IPython.core.debugger import Tracer ; debug=Tracer()

def init_particles_custom_fm(nu=0.0,mu=0.333333,a_geo  =  3.202492,b_geo  =  0.450706,mixrat_var = 'qs',numcon_var = 'qns'):
    #define new particle class by define all variables #formulated as N(D)
    #default are the variables for jplates and snowSBB size distribution
    class particle(object):
            def __init__(self, **kwargs):
                    self.__dict__.update(kwargs)
    return particle(nu_SB = nu,
                    mu_SB = mu,
                    a_geo = a_geo,
                    b_geo = b_geo,
                    mixrat_var = mixrat_var,
                    numcon_var = numcon_var)

def init_particles_custom_fd(mu=0.0,gam=1,a_ms  = 0.08 ,b_ms  =  2.22,mixrat_var = 'qs',numcon_var = 'qns'):
    #define new particle class by define all variables #formulated as N(m)
    class particle(object):
            def __init__(self, **kwargs):
                    self.__dict__.update(kwargs)
    return particle(mu = mu,
                    gam = gam,
                    a_ms = a_ms,
                    b_ms = b_ms,
                    mixrat_var = mixrat_var,
                    numcon_var = numcon_var)

def read_MC_SPlist(pathname):
    ###
    #load SP-file (keys are: "m_tot","Frim","height","d_rime","vt","xi",    "rhor","a_rime","mr_crit","diam",    "proj_A",   "mm",         "m_rime",   "m_wat")
    ###
    from netCDF4 import Dataset
    
    
    model_top=5000;heights=1;heightstep=1
    SP_file = Dataset(  pathname + '/mass2fr_' + str(tstep).zfill(4) + '-' + str(tstep_end).zfill(4) + 'min_avtstep_' + str(av_tstep) + '.ncdf',mode='r')
    #create dictionary for all variables from PAMTRA
    SP = dict()
    #if necessary change name of variables (not done yet)
    varlist = SP_file.variables
    #read Sp variables to SP dictionary
    for var in varlist:#read files and write it with different names in Data
        SP[var] = np.squeeze(SP_file.variables[var])

def calc_ND_MD(particle_class,q_h,n_tot,diam):

    '''
    calculate ND and MD for a given particle class (which must contain all paramters in N(D)=N0*D**mu*exp(-lam*D**gam) and m=a_ms*D**b_ms taking the 1st (n_tot) and 3rd (q_h) moment
    INPUT:
    particle_class: class containing relevant class parameter
    q_h: hydrometeor mass content [kg/m3]
    n_tot: hydrometeor mass concentrations [1/m3]
    diam: diameter array
    '''
    from scipy.special import gamma
    #taken from PAMTRA make_dist_params.f90
    work2 = gamma((particle_class.mu + particle_class.b_ms + 1.0) / particle_class.gam)
    work3 = gamma((particle_class.mu + 1.0) / particle_class.gam)
    lam =       (particle_class.a_ms / q_h * n_tot * work2 / work3)**(particle_class.gam / particle_class.b_ms)
    N0 = particle_class.gam * n_tot / work3 * lam**((particle_class.mu + 1.0) / particle_class.gam)
   
    N_D = N0*diam**particle_class.mu*np.exp(-lam*diam**particle_class.gam)
    M_D = N_D*particle_class.a_ms*diam**particle_class.b_ms

    return N_D,M_D

def read_and_plot_onerun(ax,pathname,filestring_MC,filestring_SB,endmin,minstep,particle_class,iwc=0,m_mean=1e10,var="Nd",unitvar="1/m3"):
    '''
    reads the McSnow (distribution...) and SB output (twomom.dat) of one run and plots them
        ax: axis to which the distributions should be plotted
        pathname: path of the directory
        filestring_MC: file with the distribution of McSnow
        filestring_SB: file with the SB moments
        minend: end of model run in minutes
        minstep: output step of model in minutes
    OPTIONAL:
        particle class: containing size dist params:
            nu,mu: parameter of the N(m) distribution N(m)=N0*m**nu*exp(-lam*m**mu)
          and m-D relations
            a_ms,b_ms fomr m(D)=a_ms*D**b_ms
        m_mean: mean mass (is passed for labelling only)
        var: which variable should be plotted; unitvar: unit of that variable
    '''

    ###plotting
    # Set the global font to be DejaVu Sans, size 10 (or any other sans-serif font of your choice!)
    rc('font',**{'family':'sans-serif','sans-serif':['DejaVu Sans'],'size':8})

    # Set the font used for MathJax - more on this later
    rc('mathtext',**{'default':'regular'})
    import matplotlib.pylab as pylab
    params = {'legend.fontsize': 6,
              'legend.handlelength': 2,
              'figure.figsize': (15, 5),
             'axes.labelsize': 'x-large',
             'axes.titlesize':'x-large',
             'xtick.labelsize':'x-large',
             'ytick.labelsize':'x-large'}

    #read the McSnow distribution file
    n_tsteps = endmin/minstep #number of timesteps in file
    MC_hei2mass = __postprocess_McSnow.read_MCdistribution(filestring_MC,n_tsteps)
   
    #convert Nd and Md from [kg/m3/ln(m)] to [kg/m3]
    MC_hei2mass["diam"] = 2.*MC_hei2mass["radius"]
    del_radius = np.diff(MC_hei2mass["radius"])
    #read_MC_SPlist(pathname)

    for key in MC_hei2mass.keys():
        if key.startswith("Nd") or key.startswith("Md"):
            MC_hei2mass[key + "_conv"] = MC_hei2mass[key][:,:-1]/del_radius*np.diff(np.log(2.*MC_hei2mass["radius"])) #/2 because of diam as x-axis
    #read 2mom-scheme bulk variables
    twomom_d = __postprocess_McSnow.read_twomom_d(filestring_SB,1) #second argument is the number of heights (should be 1)
    
    #plot distribution of different timesteps
    for i_time,time in enumerate(np.arange(0,endmin,minstep)):
        if time%20==0: #show only every ... minutes
            #get SB distributions
            twomom_d["Nd_snow"],twomom_d["Md_snow"] = __postprocess_SB.calc_distribution_from_moments_custom(twomom_d,custom_snow ,MC_hei2mass["diam"],i_time=i_time,i_height=0)
            #plot McSnow distribution
            line, = ax.loglog(MC_hei2mass["diam"][:-1],0.5*MC_hei2mass[var+ "_conv"][i_time,:],label=str(time)+"min")
            #debug()
            plot_CDF=False
            if i_time==0 and plot_CDF:#add the cdf
                ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
                #fake line for legend
                ax2.loglog(np.nan,np.nan,color=line.get_color(),label=str(time)+"min")
            if i_time==0 and plot_CDF:#add the cdf
                ax2.set_ylabel("cdf")
                ax2.loglog(MC_hei2mass["diam"][:-1],np.cumsum(MC_hei2mass[var+ "_conv"][i_time,:]*del_radius)/np.sum(MC_hei2mass[var+ "_conv"][i_time,:]*del_radius),color=line.get_color(),linewidth=0.2,label="cdf MC " + str(time)+"min")
            #plot SB distribution
            ax.loglog(MC_hei2mass["diam"][:-1],twomom_d[var + "_snow"],color=line.get_color(),linestyle='--',label="__SB"+str(time)+"min")
            if i_time==0 and plot_CDF:#add the cdf
                ax2.loglog(MC_hei2mass["diam"][:-1],np.cumsum(twomom_d[var+ "_snow"]*del_radius)/np.sum(twomom_d[var+ "_snow"]*del_radius),color=line.get_color(),linestyle='--',linewidth=0.2,label="cdf SB " + str(time)+"min")
                #debug()    
    #axis label
    ax.set_xlabel("D [m]")
    ax.set_ylabel(var + " ["+ unitvar + "]")
#return ax
    #show setted parameters
    ax.text(0.05,0.98,"mu_D={:.3f}; gam_D={:.3f};  nu_m={:.3f}; mu_m={:.3f}; \nqni={:.1e} 1/m3; iwc={:.1e}kg/m3;  m_mean={:.1e}kg  \nm={:.2f}*D**{:.2f}".format(particle_class.mu,particle_class.gam,particle_class.nu_SB,particle_class.mu_SB,iwc/m_mean,iwc,m_mean,particle_class.a_ms,particle_class.b_ms), 
            horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes)
    #debug()
    #show SB moments
    ax.text(0.05,0.1,"SB t=0: qns={:.1e}; qs={:.1e}\n    t=-1 qns={:.1e}; qs={:.1e}".format(twomom_d["qns"][0][0],twomom_d["qs"][0][0],twomom_d["qns"][-1,0],twomom_d["qs"][-1,0]), 
            horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes)
    #show MC moments
    ax.text(0.5,0.1,"MC t=0: qns={:.1e}; qs={:.1e}\n     t=-1 qns={:.1e}; qs={:.1e}".format(np.sum(MC_hei2mass["Nd_conv"][0]*del_radius),np.sum(MC_hei2mass["Md_conv"][0]*del_radius),np.sum(MC_hei2mass["Nd_conv"][-1]*del_radius),np.sum(MC_hei2mass["Md_conv"][-1]*del_radius)), 
            horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes)
    #uncomment next line if you dont use semilogy
    MC_ymax = np.nanmax(MC_hei2mass[var+ "_conv"]/2)
    SB_ymax = np.nanmax(twomom_d[var + "_snow"])
    #MC_largest_diam = MC_hei2mass["diam"][MC_hei2mass[var+ "_conv"][-1]>MC_ymax*1e-10][-1]
    #ax.set_xlim([0,1e-3]) #MC_largest_diam])
    ax.set_xlim([5e-6,0.05]) 
    if var=="Nd":
        ax.set_ylim([1e-0,1e12])
    if var=="Md":
        ax.set_ylim([1e-10,1e2])
    #ax.set_ylim([SB_ymax*1e-10,SB_ymax*2])
    #return ax,ax2
if __name__ == "__main__":
    import sys

    run_mode = int(sys.argv[1]) #0-> recompile and run 1-> run 2-> read only
    
    #define model run specifics
    endmin  = 61 #end of model run in minutes #ATTENTION:only even numbers 10,20,.. allowed; endmin is not included in output
    minstep = 1 #output step of model in minutes

    #define the parameter
    mu_arr = [5.0,1.0] #[0.,1.,2.,3.,4.,5.]
    m_mean_init_arr = [5e-11,5e-11,5e-11] #[1.19e-11,7e-11,1.19e-10]
    nrp_arr=[8388608,8388608*0.6,8388608*0.2] #in combination with m_mean_init (no own loop)
    gam=1.; 
    ########plotting 
    #define the figure grid
    plot_grid= [len(mu_arr),len(m_mean_init_arr)]
    
    #define figure for different variables
    variables=["Nd","Md"];unitvar=["1/m3","kg/m3"]
    #for var in variables:
    figNd, axesNd = plt.subplots(nrows=plot_grid[0], ncols=plot_grid[1], figsize=(plot_grid[0]*6,plot_grid[1]*4))
    figMd, axesMd = plt.subplots(nrows=plot_grid[0], ncols=plot_grid[1], figsize=(plot_grid[0]*6,plot_grid[1]*4))
    #exec("fig"+var+"=fig")
    #exec("axes"+var+"=axes")

    
    #run and plot with different settings 
    for i_mu,mu in enumerate(mu_arr):
        for i_m_mean,m_mean in enumerate(m_mean_init_arr):
            #define/read parameters describing particle properties and size distribution
            #if i_m_mean>0 or i_mu>0:
            #    continue
            nrp=nrp_arr[i_m_mean]
            custom_snow = init_particles_custom_fd(mu=mu,gam=gam,a_ms=0.08,b_ms=2.22)
            custom_snow = __postprocess_SB.convert_ND_to_Nmsingleclass(custom_snow) 
            
            #execute the run script
            if run_mode==0:
                pathname = McSnow_boxruns.run_boxi(endmin=endmin,minstep=minstep,nu=custom_snow.nu_SB,mu=custom_snow.mu_SB,iwc=m_mean*nrp,nrp=nrp,recompile=True,onlyname=False)
            elif run_mode==1:
                pathname = McSnow_boxruns.run_boxi(endmin=endmin,minstep=minstep,nu=custom_snow.nu_SB,mu=custom_snow.mu_SB,iwc=m_mean*8388608,nrp=nrp,recompile=False,onlyname=False)
            elif run_mode==2:
                pathname = McSnow_boxruns.run_boxi(endmin=endmin,minstep=minstep,nu=custom_snow.nu_SB,mu=custom_snow.mu_SB,iwc=m_mean*nrp,nrp=nrp,recompile=False,onlyname=True)
             
            #define the path of the output files
            filestring_MC=pathname + "/distribution_0000_5000.dat" 
            filestring_SB=pathname + "/twomom_d.dat" 
            
            for var,unitvar in zip(["Nd","Md"],["1/m3","kg/m3"]):
                #plot this run
                exec("axnow=axes"+ var + "[i_mu,i_m_mean]")
                read_and_plot_onerun(axnow,pathname,filestring_MC,filestring_SB,endmin,minstep,particle_class=custom_snow,m_mean=m_mean,iwc=m_mean*nrp,var=var,unitvar=unitvar)
    #for fignow,varname in zip([figNd],variables):
    for fignow,varname in zip([figNd,figMd],variables):
        plt.figure(fignow.number)
        plt.tight_layout()
        plt.plot(np.nan,np.nan,color="k",linestyle='-',label="MC")
        plt.plot(np.nan,np.nan,color="k",linestyle='--',label="SB")
        #plt.sca(ax)        
        plt.legend(loc="lower left")
        #plt.sca(ax2)        
        #plt.legend(loc="lower right")
        dir_save="/home/mkarrer/Dokumente/plots/boxi/"
        out_filestring="agg_nu_variation" + varname
        plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
        print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
        subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
