# coding: utf-8

'''
this script analyzes numerically the expectation value of E(sum(x^b)/sum(x)^b)
'''

#import packages
import numpy as np
import math
import matplotlib.pyplot as plt
import __plotting_functions
import subprocess
#from IPython.core.debugger import Tracer ; Tracer()()



#optimize the appearance of the plot (figure size, fonts)
b_array = [1.87,2.0,2.47,3.0] #tested b-parameter
number_of_plots = 2
[fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots)

#use latex
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def sum_of_xb_divided_by_b_of_sum(random_values,b):
    return sum(random_values**b)/sum(random_values)**b
    #return sum(random_values)**b
    #return sum(random_values**b)
mean=1. #mean of the distribution (actually doesnt have an effect
N_repeat_random_drawing = 200000 #number of repeating the drawing from the distributions (higher value increases the robustness of the result)
N_sum_max = 100 #maximum of evaluated N
N_array = range(1,N_sum_max+1)
N_array_numpy = np.array(N_array)

#initialize result array
#for each drawing (overwritten for different distributions)
ratio_N_drawing = np.ones([N_repeat_random_drawing,N_sum_max])*np.nan#this is the result of one random drawing
#for each distribution
ratio_monodisperse = np.ones([N_sum_max,1])*np.nan
ratio_uniform = np.ones([N_sum_max,1])*np.nan
ratio_exponential = np.ones([N_sum_max,1])*np.nan

for b in b_array:
    print "exponent ",b
    for N in range(1,N_sum_max+1):

        #monodisperse
        rand_val_monodisperse = mean*np.ones([N,1])
        ratio_monodisperse[N-1] = sum_of_xb_divided_by_b_of_sum(rand_val_monodisperse,b)

        #uniform
        #from IPython.core.debugger import Tracer ; Tracer()()
        for i_repeat in range(0,int(math.ceil(N_repeat_random_drawing/N))):
            rand_val = np.random.uniform(low=mean-0.9*mean,high=mean+0.9*mean,size=[N,1])
            ratio_N_drawing[i_repeat,N-1] = sum_of_xb_divided_by_b_of_sum(rand_val,b) 
        #mean value of the different drawings
        ratio_uniform[N-1] = np.nanmean(ratio_N_drawing[:,N-1])


        #exponential
        print  "N,N_repeat",N,int(math.ceil(N_repeat_random_drawing/N))
        for i_repeat in range(0,int(math.ceil(N_repeat_random_drawing/N))):
            rand_val= np.random.exponential(scale=mean,size=[N,1])

            ratio_N_drawing[i_repeat,N-1] = sum_of_xb_divided_by_b_of_sum(rand_val,b) 
        if N==2:
            pass #from IPython.core.debugger import Tracer ; Tracer()()
        ratio_exponential[N-1] = np.nanmean(ratio_N_drawing[:,N-1])
        
    print "monodisperse: sum(x^b)/sum(x)^b \n" , ratio_monodisperse
    print "uniform: sum(x^b)/sum(x)^b \n" , ratio_uniform
    print "exponential: sum(x^b)/sum(x)^b \n" , ratio_exponential
    print "exp/mono: ",ratio_exponential/ratio_monodisperse,
    print "exp/mono saturates at: ",np.mean(ratio_exponential[-100:-1]/ratio_monodisperse[-100:-1])
    for i_axes in range(0,2):
        if i_axes==0:
            axes[i_axes].set_xscale("linear")
            axes[i_axes].set_yscale("linear")
        elif i_axes==1:
            axes[i_axes].set_xscale("log")
            axes[i_axes].set_yscale("linear")
            
        #plot sum(x^b)/sum(x)^b for different distribution (monodisperse and exponential) and the analytical form for exponential and big N
        mono_plot = axes[i_axes].plot(N_array,ratio_monodisperse,linestyle='--') #,label="mono" + str(b))
        #plt.plot(range(1,N),ratio_uniform,color=mono_plot[0].get_color(),linestyle='-.')
        exp_plot = axes[i_axes].plot(N_array,ratio_exponential,label="b= " + str(b),color=mono_plot[0].get_color())
        
        what_function  = (1./ ( N_array_numpy**(1.-b) /ratio_exponential[:,0] ) )  #/ (1+ (N_array_numpy-1.) /   ((N_array_numpy+1.)) * (math.factorial(b)-1) )
        print "what function? ", what_function
        #from IPython.core.debugger import Tracer ; Tracer()()
        analytical_for_bigN = np.zeros_like(N_array_numpy)
        
        nominator =  np.zeros([N_sum_max]) 
        denominator = np.zeros([N_sum_max]) #np.zeros_like(N_array_numpy, np.float)
        analytical_for_bigN =  np.zeros([N_sum_max]) 
        for i_N in range(1,N_sum_max+1):
            
            nominator[i_N-1] = N_array_numpy[i_N-1]*math.gamma(1+b)
            denominator[i_N-1] = ( math.gamma(N_array_numpy[i_N-1]+b) / math.gamma(N_array_numpy[i_N-1]) )
            
            #from operator import truediv
            analytical_for_bigN[i_N-1] = math.gamma(N_array_numpy[i_N-1]+1) * math.gamma(1+b) / math.gamma(N_array_numpy[i_N-1]+b) # nominator/denominator
        print "analytical_for_bigN",analytical_for_bigN
        #from IPython.core.debugger import Tracer ; Tracer()()
        if i_axes==0:
            markevery=1
        elif i_axes==1:
            markevery=20
        #bigN_exp = axes[i_axes].plot(N_array,analytical_for_bigN,marker='x',markevery=markevery)
        bigN_exp = axes[i_axes].plot(N_array,analytical_for_bigN,marker='x',linestyle='-.',markevery=markevery,color=mono_plot[0].get_color())
        if i_axes==0:
            axes[i_axes].set_xlim([1,20])
            axes[i_axes].set_ylim([0.0,1.0])
        if i_axes==1:
            axes[i_axes].set_xlim(left=1)  
for i_axes in range(0,2):
  
    #add legend entries
    axes[i_axes].plot(np.nan,np.nan,linestyle='--',color='k',label="monodisperse: $N^(1-b)$")
    axes[i_axes].plot(np.nan,np.nan,linestyle='-',color='k',label="exponential")
    axes[i_axes].plot(np.nan,np.nan,linestyle='-.',color='k',marker='x',label="gamma(N+1)*gamma(1+b)/gamma(N+b)")
    axes[i_axes].set_ylabel(r"$E[sum(x^b)]/E[sum(x)^b]$")
    axes[i_axes].set_xlabel("N")

    axes[i_axes].legend()
plt.tight_layout()

dir_save = '/home/mkarrer/Dokumente/plots/'
out_filestring = "linear_adding_exp"
plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
plt.savefig(dir_save + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])