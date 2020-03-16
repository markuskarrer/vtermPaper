'''
this script plots the relations derived in the vterm paper and the implications for the collection kernel
'''

from IPython.core.debugger import Tracer; debug = Tracer()
import ICnew_ppropdict
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import subprocess
from scipy.special import gamma
from scipy import stats
import matplotlib.colors as colors
#get the particle properties
p = ICnew_ppropdict.init_class()

#organize the figure
fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2,ncols=2)
fig4, ax5 = plt.subplots(nrows=1,ncols=1)
#figure for bulk collision rates
fig3, ax6 = plt.subplots()
diam_lim=[1e-4,3e-2]
mass_lim=[1e-9,1e-5]
area_lim=[1e-9,1e-5] #TODO:

#initialize size array
diam = np.logspace(-5,0,200)
def_colors = plt.rcParams['axes.prop_cycle'].by_key()['color'] #default colors
for ax in [ax2,ax2,ax3,ax4,ax5,ax6]:
    ax.set_prop_cycle(color=[def_colors[1],def_colors[2],def_colors[3],def_colors[4],def_colors[5],def_colors[6],"black"])
for p_key in ["Mix2"]: #["snowSBB","aggPlate","aggNeedle","aggDendrite","aggColumn","Mix1","Mix2"]:
###plotting
    mass_p_key = p[p_key].a_ms*diam**p[p_key].b_ms #temporary: calculation of the mass array corresponding to the key
    arearatio_p_key = p[p_key].a_A*diam**p[p_key].b_A/(np.pi/4.*diam**2) #temporary: calculation of the mass array corresponding to the key
    #calcualte linear regression
    b_arat, a_arat, r_value, p_value, std_err = stats.linregress(diam,arearatio_p_key)
    #debug()
    
    


    ieffdens_p_key = 1./(mass_p_key*6./np.pi/diam**3)
    rho_l=1000.
    if p_key in ["snowSBB"]: #powerlaws
        vterm_p_key = p[p_key].a_vel*mass_p_key**p[p_key].b_vel
    else:
        vterm_p_key = p[p_key].a_Atlas-p[p_key].b_Atlas*np.exp(-p[p_key].c_Atlas*(6./np.pi*mass_p_key/rho_l)**(1./3.))
    mode_kernel= 0 #0-full kernel; 1- fix m-D (vary v-D) 2-fix v-D (vary m-D)
    if mode_kernel==0:
        pass
    elif mode_kernel==1:
        mass_p_key = p["Mix2"].a_ms*diam**p["Mix2"].b_ms #ATTENTION: test: take m-D always from Mix2
    elif mode_kernel==2:
        vterm_p_key = p["Mix2"].a_Atlas-p["Mix2"].b_Atlas*np.exp(-p["Mix2"].c_Atlas*(6./np.pi*mass_p_key/rho_l)**(1./3.))#ATTENTION: test: take v-D always from Mix2
    ax1.loglog(diam,mass_p_key,label=p_key)
    ax2.loglog(mass_p_key,diam**2)
    ax3.semilogx(diam,vterm_p_key)
    ax4.semilogx(mass_p_key,vterm_p_key)
    ax5.plot(diam,arearatio_p_key)
    ax5.plot(diam,a_arat*diam**b_arat,linestyle="--")
    calc_coll_kernel=False
    if calc_coll_kernel: 
        ###collection kernel
        #initialize collection kernel
        Kij = np.nan*np.ones((diam.shape[0],diam.shape[0]))
        for i_mass,mass_i in enumerate(mass_p_key):
            for j_mass,mass_j in enumerate(mass_p_key):
                if j_mass>i_mass:
                    continue
                Kij[i_mass,j_mass] = (diam[i_mass]+diam[j_mass])**2*(vterm_p_key[i_mass]-vterm_p_key[j_mass])

        #clear histogram
        fig2, ax5 = plt.subplots()
        bounds = np.logspace(-9,-3,20)
        cmap=plt.cm.nipy_spectral
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        #norm = colors.LogNorm(vmin=1e-9, vmax=1e-3)
        #cf = ax5.pcolormesh(diam,diam,Kij,vmin=0,vmax=1e-4,cmap="nipy_spectral_r")
        cf = ax5.pcolormesh(mass_p_key,mass_p_key,Kij,cmap=cmap,norm=norm)
        #change scale
        ax5.set_xscale('log')
        ax5.set_yscale('log')
        #limit range
        ax5.set_xlim(mass_lim) #[1e-4,5e-2])
        ax5.set_ylim(mass_lim) #[1e-4,5e-2])
        #labels
        ax5.set_xlabel("$m_{small} [kg]$")
        ax5.set_ylabel("$m_{large} [kg]$")

        #colorbar range
        cbar = fig2.colorbar(cf,ax=ax5)
        cbar.set_label("coll. kernel Kij [m3/s]")
        ax5.text(1.0, 0.0,p_key, ha='right', va='bottom', transform=ax5.transAxes,fontsize=30)
        #save figure
        if mode_kernel==0:
            fig2.savefig("/home/mkarrer/Dokumente/plots/ICnew/coll_kernel_"+ p_key  +".png")
        elif mode_kernel==1:
            fig2.savefig("/home/mkarrer/Dokumente/plots/ICnew/coll_kernel_"+ p_key  +"_vDMix2.png")
        elif mode_kernel==2:
            fig2.savefig("/home/mkarrer/Dokumente/plots/ICnew/coll_kernel_"+ p_key  +"_mDMix2.png")
    mult_by_Nsquared=False
    if mult_by_Nsquared:
        nu = 0.
        mu = 0.5
        m_mean = 1e-8 #kg
        mass_dens =1*m_mean #[g/m3] mass density
        num_dens = 1 #mass_dens/m_mean #[1/m3] number density #we get rid of this by normalizing 
        lam = (gamma((nu+1)/mu)/gamma((nu+2)/mu)*mass_dens/num_dens)**(-mu)
        N0 = mu*num_dens/gamma((nu+1)/mu)*lam**((nu+1)/mu)
        N_m = N0*mass_p_key**nu*np.exp(-lam*mass_p_key**mu)
        N2 = np.zeros((diam.shape[0],diam.shape[0]))
        KijN2 = np.zeros((diam.shape[0],diam.shape[0]))
        fig4, ax7 = plt.subplots()
        for i_mass,mass_i in enumerate(mass_p_key):
            for j_mass,mass_j in enumerate(mass_p_key):
                if j_mass>i_mass:
                    continue
                
        
                N2[i_mass,j_mass] += N_m[i_mass]*N_m[j_mass]
        

        for i_mass,mass_i in enumerate(mass_p_key):
            for j_mass,mass_j in enumerate(mass_p_key):
                #KijN2[i_mass,j_mass] = Kij[i_mass,j_mass]*N2[i_mass,j_mass] 
                KijN2[i_mass,j_mass] = N2[i_mass,j_mass] 
        #debug()
        cfN2 = ax7.contour(mass_p_key,mass_p_key,KijN2) #,cmap=cmap,norm=norm)
        #change scale
        ax7.set_xscale('log')
        ax7.set_yscale('log')
        #limit range
        ax7.set_xlim(mass_lim) #[1e-4,5e-2])
        ax7.set_ylim(mass_lim) #[1e-4,5e-2])
        #labels
        ax7.set_xlabel("$m_{small} [kg]$")
        ax7.set_ylabel("$m_{large} [kg]$")

        #colorbar range
        cbar = fig4.colorbar(cf,ax=ax7)
        cbar.set_label("coll. kernel Kij*N**2 [s/m3]")
        ax7.text(1.0, 0.0,p_key, ha='right', va='bottom', transform=ax7.transAxes,fontsize=30)
        #save figure
        if mode_kernel==0:
            fig4.savefig("/home/mkarrer/Dokumente/plots/ICnew/coll_kernel_N2_"+ p_key  +".png")
        elif mode_kernel==1:
            fig4.savefig("/home/mkarrer/Dokumente/plots/ICnew/coll_kernel_N2_"+ p_key  +"_vDMix2.png")
        elif mode_kernel==2:
            fig4.savefig("/home/mkarrer/Dokumente/plots/ICnew/coll_kernel_N2_"+ p_key  +"_mDMix2.png")
         


    calc_bulk_coll_rates=False
    if calc_bulk_coll_rates:
        m_array_bulk=mass_p_key[::10]
        bulk_rates = np.ones_like(m_array_bulk)
        for i_m_array_bulk,m_mean in enumerate(m_array_bulk): #loop over differnt mean mass
            nu = 0.
            mu = 0.5
            mass_dens =1*m_mean #[g/m3] mass density
            num_dens = 1 #mass_dens/m_mean #[1/m3] number density #we get rid of this by normalizing 
            lam = (gamma((nu+1)/mu)/gamma((nu+2)/mu)*mass_dens/num_dens)**(-mu)
            N0 = mu*num_dens/gamma((nu+1)/mu)*lam**((nu+1)/mu)
            
            N_m = N0*mass_p_key**nu*np.exp(-lam*mass_p_key**mu)
            bulk_rates[i_m_array_bulk] = 0 
            for i_mass,mass_i in enumerate(mass_p_key):
                for j_mass,mass_j in enumerate(mass_p_key):
                    if j_mass>i_mass:
                        continue
                    if Kij[i_mass,j_mass]>0:
                        bulk_rates[i_m_array_bulk] += N_m[i_mass]*N_m[j_mass]*Kij[i_mass,j_mass]
                    if i_mass==10 & j_mass==10:
                        print "m_mean,i_mass,j_mass,mass_i,mass_j,N_m[i_mass],N_m[j_mass],Kij[i_mass,j_mass]"
                        print m_mean,i_mass,j_mass,mass_i,mass_j,N_m[i_mass],N_m[j_mass],Kij[i_mass,j_mass]
        
        ax6.loglog(m_array_bulk[:-1],bulk_rates[:-1],label=p_key)#/np.diff(m_array_bulk))

        ax6.set_xlabel("mean mass m$_{mean}$ [kg]")
        ax6.set_ylabel("normalized bulk coll. rate [m3/s/kg]")
        ax6.set_xlim([1e-12,mass_lim[1]]) #[1e-4,5e-2])
if calc_bulk_coll_rates:
    ax6.legend()
    #save figure 
    if mode_kernel==0:
        fig3.savefig("/home/mkarrer/Dokumente/plots/ICnew/bulk_collision_rates_.png")
    elif mode_kernel==1:
        fig3.savefig("/home/mkarrer/Dokumente/plots/ICnew/bulk_collision_rates_vDMix2.png")
    elif mode_kernel==2:
        fig3.savefig("/home/mkarrer/Dokumente/plots/ICnew/bulk_collision_rates_mDMix2.png")
##labelling
#ax1.set_xlabel("Diameter [m]")
#ax2.set_xlabel("mass [kg]")
ax3.set_xlabel("$D_{max}$ [m]")
ax4.set_xlabel("mass m [kg]")
ax5.set_xlabel("$D_{max}$ [m]")
#ylabel
ax1.set_ylabel("mass m [kg]")
ax2.set_ylabel("$D_{max}^2$ [m2]")
ax3.set_ylabel("v$_{term}$ [m/s]")
ax5.set_ylabel("area ratio $A_r$")

##limits
#xlims
ax1.set_xlim(diam_lim)
ax2.set_xlim(mass_lim)
ax3.set_xlim(diam_lim)
ax4.set_xlim(mass_lim)
ax5.set_xlim([1e-4,1e-2]) #diam_lim)
#ylims
ax1.set_ylim(mass_lim)
ax2.set_ylim([1e-8,1e-3])
ax3.set_ylim([0,2.5])
ax4.set_ylim([0,2.5])
ax5.set_ylim([0,0.8])
ax5.set_yticks([0.2,0.4,0.6,0.8])



ax1.legend()

plt.tight_layout()
fig1.savefig('/home/mkarrer/Dokumente/plots/ICnew/pprop_vtermpaper.pdf', dpi=301)

fig4.savefig('/home/mkarrer/Dokumente/plots/ICnew/asratio.pdf')

#if mode_kernel==0:
#    subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/ICnew/pprop_vtermpaper.pdf'])
