
import subprocess
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import __postprocess_McSnow
import __postprocess_SB
import oned_SB_tests_run
import vtermstudy_coeffs
from matplotlib import rc 
'''
this code executes McSnow_boxruns.py with different settings (particle properties,size distribution parameter, initialization (nrp0,iwc), supersaturation) and visualizes its output
'''
from IPython.core.debugger import Tracer ; debug=Tracer()

def read_and_plot_onerun(ax,pathname_dic,endmin,minstep,particle_class,iwc=0,m_mean=1e10,nrp=8388608,monotype_list=["plate"],first_plot=False):
    '''
    reads the McSnow (distribution...) and SB output (twomom.dat) of one run and plots them
        ax: axis to which the distributions should be plotted
        pathname_dic: contains path of the directories
        filestring_SB: file with the SB moments
        minend: end of model run in minutes
        minstep: output step of model in minutes
    OPTIONAL:
        particle class: containing size dist params:
            nu,mu: parameter of the N(m) distribution N(m)=N0*m**nu*exp(-lam*m**mu)
          and m-D relations
            a_ms,b_ms fomr m(D)=a_ms*D**b_ms
        m_mean: mean mass (is passed for labelling only)
        nrp: number concentration
        first_plot: put special info (legend,..) on first plot
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

   
    axnum = ax.twiny() #add an axes for the number (flux)
    for monotype in monotype_list:
        #read 2mom-scheme bulk variables
        #debug()
        twomom_d = __postprocess_McSnow.read_twomom_d(pathname_dic[monotype] + "/twomom_d.dat" 
     ,250) #second argument is the number of heights (should be 1)
        #add height
        twomom_d["height"]= 5000-np.cumsum(twomom_d["dz"][0])
        axqs, = ax.semilogx(twomom_d["qs"][-1],twomom_d["height"]) #twomom_d["qns"][-1] -> always take last timestep (in equilibrium)
        axnum.semilogx(twomom_d["qns"][-1],twomom_d["height"],linestyle='--',color=axqs.get_color()) #twomom_d["qns"][-1] -> always take last timestep (in equilibrium)
        if first_plot:
            ax.plot(np.nan,np.nan,color=axqs.get_color(),label=monotype)
    if first_plot:
        ax.plot(np.nan,np.nan,color="k",linestyle='-',label="qs")
        ax.plot(np.nan,np.nan,color="k",linestyle='--',label="qns")
        ax.legend(loc="lower left")
    #axes labels
    ax.set_ylabel("height [m]")
    ax.set_xlabel("qi [kg/kg]")
    axnum.set_xlabel("qni [1/kg]")
    #debug() 
    #show setted parameters
    #show SB moments #TODO: show SB surface vars?
    #ax.text(0.05,0.1,"SB t=0: qns={:.1e}; qs={:.1e}\n    t=-1 qns={:.1e}; qs={:.1e}".format(twomom_d["qns"][0][0],twomom_d["qs"][0][0],twomom_d["qns"][-1,0],twomom_d["qs"][-1,0]), 
    #        horizontalalignment='left',
    #        verticalalignment='top', transform=ax.transAxes)
    #TODO: set limits
    #ax.set_xlim([5e-6,0.05]) 
    #if var=="Nd":
    #    ax.set_ylim([1e-0,1e12])
    #if var=="Md":
    #    ax.set_ylim([1e-10,1e2])
    #ax.set_ylim([SB_ymax*1e-10,SB_ymax*2])

if __name__ == "__main__":
    import sys

    run_mode = int(sys.argv[1]) #0-> recompile and run 1-> run 2-> read only
    
    #define model run specifics
    endmin  = 600 #end of model run in minutes #ATTENTION:only even numbers 10,20,.. allowed; endmin is not included in output
    minstep = 1 #output step of model in minutes
    

    #define the parameter
    mu_arr = [1.0,0.0] #[0.,1.,2.,3.,4.,5.]
    m_mean_init_arr = [5e-11,5e-10] #[1.19e-11,7e-11,1.19e-10]
    nrp_arr=[1e6,1e5] #in combination with m_mean_init (no own loop)
    gam=1.; 
    ########plotting 
    #define the figure grid
    plot_grid= [len(mu_arr),len(m_mean_init_arr)]
    
    #define figure for different variables
    #for var in variables:
    fig, axes = plt.subplots(nrows=plot_grid[0], ncols=plot_grid[1], figsize=(plot_grid[0]*6,plot_grid[1]*4))

    
    #run and plot with different settings 
    for i_mu,mu in enumerate(mu_arr):
        for i_m_mean,m_mean in enumerate(m_mean_init_arr):
            #define/read parameters describing particle properties and size distribution
            if i_m_mean>0 or i_mu>0:
                continue
                first_plot=False
            else:
                first_plot=True
            nrp=nrp_arr[i_m_mean]
            #custom_snow = vtermstudy_coefqns.return_class(monotype=monotype,mu=mu,gam=gam)   
            monotype_list=["plate","mixcolumndend","powerlawOLD"]
            pathname_dic=dict() #save the folder names of the runs
            for monotype in monotype_list:
                if monotype=="powerlawOLD":
                    custom_snow = vtermstudy_coeffs.return_class(monotype=monotype,nu_SB=0.0,mu_SB=0.0) 
                else:
                    custom_snow = vtermstudy_coeffs.return_class(monotype=monotype,mu=mu,gam=gam) 

                #execute the run script
                if run_mode==0:
                    pathname = oned_SB_tests_run.run_1d(monomer_type=monotype,endmin=endmin,minstep=minstep,nu=custom_snow.nu_SB,mu=custom_snow.mu_SB,iwc=m_mean*nrp,nrp=nrp,recompile=True,onlyname=False)
                elif run_mode==1:
                    print "do not use this because recompiling is esential for modifications in the SB scheme particle properties for this script"; sys.exit(1)
                    pathname = oned_SB_tests_run.run_1d(monomer_type=monotype,endmin=endmin,minstep=minstep,nu=custom_snow.nu_SB,mu=custom_snow.mu_SB,iwc=m_mean*nrp,nrp=nrp,recompile=False,onlyname=False)
                elif run_mode==2:
                    pathname = oned_SB_tests_run.run_1d(monomer_type=monotype,endmin=endmin,minstep=minstep,nu=custom_snow.nu_SB,mu=custom_snow.mu_SB,iwc=m_mean*nrp,nrp=nrp,recompile=False,onlyname=True)
            
                #define the path of the output files
                #pathname_splitted=pathname.split('_')
                #pathname="_".join(pathname_splitted[0:-1])
                pathname_dic[monotype]=pathname


             
            read_and_plot_onerun(axes[i_mu,i_m_mean],pathname_dic,endmin,minstep,particle_class=custom_snow,m_mean=m_mean,nrp=nrp,iwc=m_mean*nrp,monotype_list=monotype_list,first_plot=True)
    
    #add additional labels, save and open
    plt.tight_layout()
    dir_save="/home/mkarrer/Dokumente/plots/1d_SB_test/"
    out_filestring="profiles"
    plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
    print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
    subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
