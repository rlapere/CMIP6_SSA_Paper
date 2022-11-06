import pymannkendall as mk
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,BoundaryNorm
import netCDF4
import numpy as np
import proplot as plot
import pandas as pd
import cartopy.crs as ccrs
import cartopy
import matplotlib as mpl
from scipy.signal import savgol_filter
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
from matplotlib.collections import PathCollection
from matplotlib.path import Path

si0_ = netCDF4.Dataset('trend_siconc.nc','r')
lonssi = si0_.variables['lon'][:]
latssi = si0_.variables['lat'][:]
trendssi = si0_.variables['trend'][:,:]*12*10
vmani = np.nanpercentile(np.abs(trendssi[latssi>45].flatten()),99)
vmasi = np.nanpercentile(np.abs(trendssi[latssi<-50].flatten()),99)

w0_ = netCDF4.Dataset('trend_sfcwind.nc','r')
lonw = w0_.variables['lon'][:]
latw = w0_.variables['lat'][:]
trendw = w0_.variables['trend'][:,:]*12*10
vmaw = np.nanpercentile(np.abs(trendw[latw<-50].flatten()),99)

fig, axs = plot.subplots(ncols=2,nrows=2,proj={1:'npaeqd',2:'spaeqd',3:'npaeqd',4:'spaeqd'},sharey=False,sharex=False)
axs.format(suptitle='Trends in sea-ice and surface wind speed 1951$-$2014',coast=True, latlines=10,coastcolor='black',land=True,landcolor='white')

a=axs[0].pcolormesh(lonssi,latssi,trendssi,cmap='ColdHot',vmin=-vmasi,vmax=vmasi,extend='both')
axs[0].format(boundinglat=60,title='Arctic sea-ice')

b=axs[1].pcolormesh(lonssi,latssi,trendssi,cmap='ColdHot',vmin=-vmasi,vmax=vmasi,extend='both')
axs[1].format(boundinglat=-60,title='Antarctic sea-ice')
axs[1].colorbar(b,shrink=0.8,label='%$\,$decade$^{-1}$')

axs[2].pcolormesh(lonw,latw,trendw,cmap='ColdHot',vmin=-vmaw,vmax=vmaw,extend='both')
axs[2].format(boundinglat=60,title='Arctic wind speed')

c=axs[3].pcolormesh(lonw,latw,trendw,cmap='ColdHot',vmin=-vmaw,vmax=vmaw,extend='both')
axs[3].format(boundinglat=-60,title='Antarctic wind speed')
axs[3].colorbar(c,shrink=0.8,label='m$\,$s$^{-1}$$\,$decade$^{-1}$')

plt.show()

si0_ = netCDF4.Dataset('trend_mmrss_rel_ssp585_siconc.nc','r')
lonssi = si0_.variables['lon'][:]
latssi = si0_.variables['lat'][:]
trendssi = si0_.variables['trend'][:,:]*12*10
vmani = np.nanpercentile(np.abs(trendssi[latssi>45].flatten()),99)
vmasi = np.nanpercentile(np.abs(trendssi[latssi<-50].flatten()),99)
vmasi = np.nanmax(np.abs(trendssi.flatten()))


fig, axs = plot.subplots(ncols=1,nrows=1,proj={1:'robin'},sharey=False,sharex=False)
axs.format(suptitle='Trends in sea-ice in SSP585 2015$-$2100',coast=True, latlines=10,coastcolor='black',land=True,landcolor='white')

a=axs[0].pcolormesh(lonssi,latssi,trendssi,cmap='Blues6',vmin=-vmasi,vmax=0,extend='both')


fig, axs = plot.subplots(ncols=1,nrows=2,proj={1:'npaeqd',2:'spaeqd'},sharey=False,sharex=False)
axs.format(suptitle='Trends in sea-ice in SSP585 2015$-$2100',coast=True, latlines=10,coastcolor='black',land=True,landcolor='white')

a=axs[0].pcolormesh(lonssi,latssi,trendssi,cmap='ColdHot',vmin=-vmasi,vmax=vmasi,extend='both')
axs[0].format(boundinglat=50,title='Arctic sea-ice')

b=axs[1].pcolormesh(lonssi,latssi,trendssi,cmap='ColdHot',vmin=-vmasi,vmax=vmasi,extend='both')
axs[1].format(boundinglat=-50,title='Antarctic sea-ice')
axs[1].colorbar(b,shrink=0.8,label='%$\,$decade$^{-1}$')

plt.show()
