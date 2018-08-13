#import packages
import numpy as np
import matplotlib.pyplot as plt

def colbar_lowest_white(cmap,lenbounds=100):
    #adds a white color for the lowest values to the predefined colorbar 
    exec("colors1=plt.cm." + cmap + "(np.linspace(0.0,1.0,lenbounds))")
    colors2=plt.cm.Greys(0)
    colors_merged=np.vstack((colors2,colors1))
    new_cmap = colors.ListedColormap([colors_merged[i] for i in range(0,lenbounds)])
    return new_cmap