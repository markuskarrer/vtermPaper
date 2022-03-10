'''
plot the precipitation field from ICON output
'''

import numpy as np
from matplotlib import colors, cm, pyplot as plt
from matplotlib.ticker import NullFormatter, MaxNLocator
from netCDF4 import Dataset		#reading ncdf
import matplotlib.colors
import datetime #for conversion between unixtime and calender time
import time
from IPython.core.debugger import Tracer ; debug = Tracer()
from subprocess import call
import os

class Domain():

    def __init__(self,name,center_lon,center_lat,number_lon,number_lat,domain_rad):
    
        self.name = name
        self.center_lon = center_lon # Longitude of domain center
        self.center_lat = center_lat # Latitude of domain center
        self.number_lon = number_lon # Number of points in the longitude
        self.number_lat = number_lat # Number of points in the latitude
        self.domain_rad = domain_rad # Radius of the domain in degrees (about an isocircle)

    def getLonFirst(self):
        return self.center_lon - self.domain_rad

    def getLatFirst(self):
        return self.center_lat - self.domain_rad

    def getLonIncrement(self):
        return 2.*self.domain_rad/self.number_lon

    def getLatIncrement(self):
        return 2.*self.domain_rad/self.number_lat

    def getLonNumber(self):
        return self.number_lon
    
    def getLatNumber(self):
        return self.number_lat

def writeGridDef(dom):

    fname = "latlongriddef_"+dom.name

    with open(fname,"w") as tf:
        tf.write("gridtype = lonlat \n")
        tf.write("xsize = {0} \n".format(dom.getLonNumber()))
        tf.write("ysize = {0} \n".format(dom.getLatNumber()))
        tf.write("xfirst = {0} \n".format(dom.getLonFirst()))
        tf.write("xinc = {0} \n".format(dom.getLonIncrement()))
        tf.write("yfirst = {0} \n".format(dom.getLatFirst()))
        tf.write("yinc = {0} \n".format(dom.getLatIncrement()))
    
    return fname


import argparse
parser =  argparse.ArgumentParser(description='plot hydrometeors')
parser.add_argument('--date', nargs=1, help='gimme datestring in the format YYYYMMDD')
parser.add_argument('-int','--interpolate', nargs=1, help='interpolate latlon')
parser.add_argument('-exp','--experiment', nargs=1,
                         help='gimme the name of the sensitivity experiment')
args = parser.parse_args()
experiment = args.experiment[0]
date = args.date[0]
if args.interpolate[0]=="0":
    latlon_interp = False
else:
    latlon_interp = True

savepath = "/home/mkarrer/Dokumente/plots/ICON/precip_field/" + date + "/"

if not os.path.exists(savepath):
    os.makedirs(savepath)

if latlon_interp:
    DOM1 = Domain("DOM1",6.4145,50.9089,200,200,0.75)

    filedom1 = writeGridDef(DOM1)

    runFld = '/data/optimice/ICON_output/' + date + '_110km_' + experiment+ '/' 

    file = runFld + "626_10m_precip_DOM01_ML_" + date + "T000000Z"

    res = call(["cdo","-P","2","-r","remapnn,"+filedom1,"-selvar,rain_gsp_rate,snow_gsp_rate,ice_gsp_rate,graupel_gsp_rate,hail_gsp_rate",file + '.nc',file + '_latlon.nc'])

input_folder='/data/optimice/ICON_output/' + date + '_110km_' + experiment + '/'


input_file='626_10m_precip_DOM01_ML_' + date + 'T000000Z_latlon.nc'
filepath = input_folder + input_file

#load netcdf-file
out=Dataset(filepath)
print filepath

for i_timestep,timestep in enumerate(out["time"][:]): #timestep in file

    #debug()

    #var_range_y=[90,410] #limit variable in y-dimension because the domain is smaller than the interpolation
    #rain=out.variables['rain_gsp'][i_timestep,var_range_y[0]:var_range_y[1],:]
    rain=out.variables['rain_gsp_rate'][i_timestep,:,:]
    ice=out.variables['ice_gsp_rate'][i_timestep,:,:]
    snow=out.variables['snow_gsp_rate'][i_timestep,:,:]
    graupel=out.variables['graupel_gsp_rate'][i_timestep,:,:]
    hail=out.variables['hail_gsp_rate'][i_timestep,:,:]

    tot_prec=rain+ice+snow+graupel+hail

    #convert from kg m-3 s-1 to kg m-3 h-1
    tot_prec=tot_prec*3600.

    # Set up the size of the figure
    #	fig.clf()
    fig = plt.figure(1, figsize=(11,9))
    ax = plt.axes()

    ##plot-properties
    mincol=0;maxcol=2*10**-3#set max and min of colors
    colticksmin=0;colticksmax=maxcol*1.001;colticks=maxcol/10#specify colorbar ticks
    colorbar_label='m-3'#set colorbar label
    #colorbar from: http://jjhelmus.github.io/blog/2013/09/17/plotting-nsw-precipitation-data/
    nws_precip_colors = [
        "#d3d3d3", #"#fdfdfd",   # 10.00+
        "#04e9e7",  # 0.01 - 0.10 inches
        "#019ff4",  # 0.10 - 0.25 inches
        "#0300f4",  # 0.25 - 0.50 inches
        "#02fd02",  # 0.50 - 0.75 inches
        "#01c501",  # 0.75 - 1.00 inches
        "#008e00",  # 1.00 - 1.50 inches
        "#fdf802",  # 1.50 - 2.00 inches
        "#e5bc00",  # 2.00 - 2.50 inches
        "#fd9500",  # 2.50 - 3.00 inches
        "#fd0000",  # 3.00 - 4.00 inches
        "#d40000",  # 4.00 - 5.00 inches
        "#bc0000",  # 5.00 - 6.00 inches
        "#f800fd",  # 6.00 - 8.00 inches
        "#9854c6",  # 8.00 - 10.00 inches
        "#fdfdfd"   # 10.00+
    ]
    precip_colormap = matplotlib.colors.ListedColormap(nws_precip_colors)
    levels = [0, 0.1, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0,
              6.0, 8.0, 10.,20.]
    #levels = [0, 0.001, 0.01, 0.1, 0.25 ,0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0,
    #	  6.0, 8.0, 10.,20.]
    cmap_name='Blues' #'GnBu' #gist_earth_r' #'YlGnBu'
    #modify colormap
    norm = matplotlib.colors.BoundaryNorm(levels, 16)
    cmap 	= 	plt.get_cmap(cmap_name)
    ###plot command
    #[lon_mesh,lat_mesh] = np.meshgrid(out.variables['lon'],out.variables['lat'][var_range_y[0]:var_range_y[1]])
    [lon_mesh,lat_mesh] = np.meshgrid(out.variables['lon'],out.variables['lat'])
    pcol	=	ax.pcolormesh(lon_mesh,lat_mesh,tot_prec,rasterized=True,norm=norm,cmap=precip_colormap)
    plt.xlabel('Longitude / $^\circ$ E')
    plt.ylabel('Latitude / $^\circ$ N')
    #plt.xlim([5.23,7.58])
    ax.set_xlim([5.70,7.2])
    plt.ylim([50.4,51.4])
    #set up colorbar
    col = fig.colorbar(pcol,ticks=levels,format='%5.2f') #run output/scripts/colbar_for_prec_field.py to get horizontal colorbar
    col.ax.set_ylabel('Precipitation rate / mm h-1', rotation=90)#colorbar label
    #col	=	plt.colorbar(pcol,ticks=np.arange(colticksmin,colticksmax,colticks))
    #add marker for JOYCE
    plt.scatter(6.4145,50.9089,s=200,marker='x',color='black')
    ax.set_title(experiment)
    #print text to know on which run/timestep you are referring
    time_HHMM = time.strftime('%H:%M',time.gmtime(out["time"][i_timestep].data*86400.)) + ' UTC'
    print time_HHMM
    plt.text(7.25,50.32, time_HHMM,{'color': 'k', 'fontsize': 30})
    filename='prec_field' + input_folder.translate(None, '/') + input_file.translate(None, '/.') + '_' + time.strftime('%H:%M',time.gmtime(out["time"][i_timestep].data*86400.)) + 'UTC'
    # Save to a File
    #make labels etc larger
    matplotlib.rcParams.update({'font.size': 25})


    #fig.savefig(savepath + filename + '.pdf',format = 'pdf')
    fig.savefig(savepath + filename + '.png',format = 'png')
    fig.clf()#avoid overlaying of previous timesteps
    print 'file is at: ' + savepath + filename + '.pdf'
