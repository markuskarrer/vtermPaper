'''
plot the precipitation field from RADOLAN
RUN WITH wradlib (conda activate wradlib)

download files with e.g. wget https://opendata.dwd.de/climate_environment/CDC/grids_germany/5_minutes/radolan/reproc/2017_002/bin/2018/YW2017.002_201811.tar -P 2018
extract twice with e.g tar -xf <files>.tar -C <make this dir first>
'''

import wradlib as wrl
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from IPython.core.debugger import Tracer ; debug = Tracer()
import os

#get arguments
import argparse
parser =  argparse.ArgumentParser(description='plot RADOLAN data')
parser.add_argument('-d','--date', nargs=1, help='gimme datestring in the format YYYYMMDD')
args = parser.parse_args()
date = args.date[0]

#get all infos from date for the RADOLAN folder structure
year = date[0:4]
year_short = year[2:4]
month = date[4:6]
day = date[6:8]

# Get coordinates
radolan_grid_xy = wrl.georef.get_radolan_grid(900,900)
radolan_egrid_xy = wrl.georef.get_radolan_grid(1500,1400)
radolan_wgrid_xy = wrl.georef.get_radolan_grid(1100, 900)
x = radolan_grid_xy[:,:,0]
y = radolan_grid_xy[:,:,1]

xe = radolan_egrid_xy[:,:,0]
ye = radolan_egrid_xy[:,:,1]

xw = radolan_wgrid_xy[:,:,0]
yw = radolan_wgrid_xy[:,:,1]

#radolan_grid = wrl.georef.get_radolan_grid()
radolan_grid_coord = wrl.georef.get_radolan_grid(1100, 900, wgs84=True)
longitude = radolan_grid_coord[:,:,0]
latitude = radolan_grid_coord[:,:,1]

def read_radolan(radfile):
    return wrl.io.read_radolan_composite(radfile)

def plot_radolan(ax,data, attrs, grid, clabel=None):

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
    precip_colormap = mpl.colors.ListedColormap(nws_precip_colors)
    levels = [0, 0.1, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0,
              6.0, 8.0, 10.,20.]
    #levels = [0, 0.001, 0.01, 0.1, 0.25 ,0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0,
    #	  6.0, 8.0, 10.,20.]
    #modify colormap
    norm = mpl.colors.BoundaryNorm(levels, 16)

    x = grid[:,:,0]
    y = grid[:,:,1]
    pm = ax.pcolormesh(x, y, data*10., norm=norm,cmap=precip_colormap)
    cb = fig.colorbar(pm,ticks=levels)
    cb.set_label(clabel)
    ax.set_xlabel('Longitude / $^\circ$ E')
    ax.set_ylabel('Latitude / $^\circ$ N')
    #ax.set_title('{0} Product\n{1}'.format(attrs['producttype'],
    #                                   attrs['datetime'].isoformat()))
    ax.set_title('{0} Product'.format(attrs['producttype']))
    #plt.xlim((x[0,0],x[-1,-1]))
    #plt.ylim((y[0,0],y[-1,-1]))
    ax.set_xlim([5.70,7.2])
    ax.set_ylim([50.4,51.4])
    #ax.grid(color='r')

    return ax

#path for output
savepath = "/home/mkarrer/Dokumente/plots/ICON/precip_field_radolan/" + date + "/"

#create path if not there yet
if not os.path.exists(savepath):
    os.makedirs(savepath)

for hour in range(0,24):
    for minute in range(0,56,5):
        time_HHMM='{:02d}'.format(hour) + '{:02d}'.format(minute)
        time_HHMMnice='{:02d}'.format(hour) + ':' + '{:02d}'.format(minute)

        fig = plt.figure(1, figsize=(11,9))
        ax = plt.axes()

        radolanPath = '/data/optimice/ICON_output/radolan/' + year + '/' + year + month + '/' + day + '/'
        radolanFile = "raa01-yw2017.002_10000-" + year_short + month + day + time_HHMM + "-dwd---bin"

        print(radolanFile)

        #read file
        data, attrs = read_radolan(radolanPath + radolanFile)

        #plot
        plot_radolan(ax,data, attrs, radolan_grid_coord, clabel='Precipitation rate / mm h-1')
        
        #add marker for JOYCE
        ax.scatter(6.4145,50.9089,s=200,marker='x',color='black')
        #add timestamp
        plt.text(7.25,50.32, time_HHMMnice + ' UTC',{'color': 'k', 'fontsize': 30})

        mpl.rcParams.update({'font.size': 25})
        fig.savefig(savepath + time_HHMM + '.png',format = 'png')
        fig.clf()
