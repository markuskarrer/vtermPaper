
import subprocess
import os
import sys
'''
this code runs the box-model of McSnow with different settings (particle properties,size distribution parameter, initialization (nrp0,iwc), supersaturation) 
'''
from IPython.core.debugger import Tracer ; debug=Tracer()

#general settings
MC_dir="/home/mkarrer/Dokumente/McSnowboxtests/mcsnow/"
ORG="/home/mkarrer/Dokumente/bash_organizing_code"

def run_1d(recompile=True,monomer_type="plate",endmin=60,minstep=5,nu=0.0,mu=1./3,iwc=0.002,nrp=8388608,onlyname=False):
    '''
    monomer_type: e.g. plate, ... #ATTENTION: for using other monotypes vterm Atlas have to be defined in mo_modify_MC.src and a new type has to be defined in mo_2mom_mcrph_main.f90
    minend: end of model run in minutes
    minstep: output step of model in minutes
    #distribution parameters in N(m)=N0*m**nu*(-lam*m**mu)
    nu: width parameter
    mu: exponential parameter

    onlyname: do not run McSnow, but just retrieve the path with this settings
    '''
    
    #convert nu and mu to strings
    nu_str="{:.6f}".format(nu)
    mu_str="{:.6f}".format(mu)

    ####execute bash commands
    #modify the src files to use desired particle properties
    os.chdir(ORG)
    

    if monomer_type in ["plate","dendrite","column","needle","rosette","mixcolumndend"]:
        call_out = subprocess.call("bash modify_SB_src.sh" + " nu" + nu_str+ "_mu" + mu_str + "_AtlasJ" + monomer_type+ " " + "snowinit "+ "box",shell=True); 
    else:
        call_out = subprocess.call("bash modify_SB_src.sh" + " nu" + nu_str+ "_mu" + mu_str + "_" + monomer_type+ " " + "snowinit "+ "box",shell=True); #ATTENTION the nu/mu-modification only have an effect for the snowj..nonsphere classes 

    if call_out!=0: sys.exit(1)

    #recompile the code
    if recompile:
        os.chdir(MC_dir)
        call_out = subprocess.call("make clean".split())
        if call_out!=0: sys.exit(1)

        call_out = subprocess.call("make release_omp".split())
        if call_out!=0: sys.exit(1)

    os.chdir(MC_dir + "run/")
    #get the output string from the box_model
    process= subprocess.Popen(("bash 1d_SDA_markus onlyname " + str(endmin) + " " + str(minstep) + " "+ nu_str +" "+ mu_str + " " + str(iwc) + " " + str(int(nrp))+ " " + monomer_type).split(),stdout=subprocess.PIPE) 
    out, err = process.communicate()
    
    pathname = out[:-1]
    if onlyname: #exit here if only the path is requested
        return pathname
        sys.exit(0)
    #run the box model
    subprocess.call("bash 1d_SDA_markus run " + str(endmin) + " " + str(minstep) +" "+  nu_str +" "+ mu_str + " " + str(iwc) + " " + str(int(nrp))+ " " + monomer_type ,shell=True) 

    return pathname 

if __name__ == "__main__":
    
    #first argument in calling script determines if it should be recompiled
    recompile=(sys.argv[1]=='1')
    #debug() 
    run_boxi(recompile=recompile)
    run_boxi(nu=5)
    run_boxi(nu=10)
