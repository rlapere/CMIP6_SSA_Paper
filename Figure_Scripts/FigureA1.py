########################
# Remy Lapere (03/11/2022)
# Figure A1 for Lapere et al., 2022 - CMIP6 SSaer
# input files for the script are from the MACv2 climatology (Kinne, 2019)
########################

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import proplot as plot
import pandas as pd
from matplotlib.colors import LogNorm
from scipy.stats import pearsonr
import matplotlib as mpl

fld = '/Users/rlapere/Desktop/CMIP_SSA/gt_aodSS_550.nc'
md = netCDF4.Dataset(fld,'r')
modis = md.variables['com_AOD'][:]
modisavg = md.variables['com_AOD_ann'][:]
latsm = md.variables['lat'][:]
lonsm = md.variables['lon'][:]
ddts = md.variables['time'][:]

fld1 = '/Users/rlapere/Desktop/CMIP_SSA/gt_aodDU_550.nc'
md1 = netCDF4.Dataset(fld1,'r')
modis1 = md1.variables['com_AOD'][:]
modisavg1 = md1.variables['com_AOD_ann'][:]
latsm1 = md1.variables['lat'][:]
lonsm1 = md1.variables['lon'][:]
ddts1 = md1.variables['time'][:]

frac = '/Users/rlapere/Desktop/CMIP_SSA/gt_ff_0550nm.nc'
ff_ = netCDF4.Dataset(frac,'r')
ff = ff_.variables['aer_data'][:]*100
latsf = ff_.variables['lat'][:]
lonsf = ff_.variables['lon'][:]

fld11 = '/Users/rlapere/Desktop/CMIP_SSA/gt_aodTO_550.nc'
md11 = netCDF4.Dataset(fld11,'r')
modis11 = md11.variables['com_AOD_ann'][:]
modrat = np.mean(modis/(modis+modis1),axis=0)*100

vmi,vma=0,100

f, axs = plot.subplots(ncols=3,nrows=2,proj={1:'npaeqd',2:'npaeqd',3:'npaeqd',4:'spaeqd',5:'spaeqd',6:'spaeqd'})
axs.format(coast=True, latlines=10,coastcolor='black',land=True,landcolor='white',suptitle='AOD climatology from MACv2')

m=axs[0].pcolormesh(lonsm,latsm,modrat,cmap='ColdHot',vmin=vmi,vmax=vma)
m1=axs[3].pcolormesh(lonsm,latsm,modrat,cmap='ColdHot',vmin=vmi,vmax=vma)
axs[0].format(boundinglat=50)
axs[3].format(boundinglat=-50)
cbar = axs[0].colorbar(m,label='%')
axs[0].set_title('Sea-salt AOD/(sea-salt+dust AOD)')


m=axs[1].pcolormesh(lonsf,latsf,np.mean(ff,axis=0),cmap='ColdHot',vmin=vmi,vmax=100)
m1=axs[4].pcolormesh(lonsf,latsf,np.mean(ff,axis=0),cmap='ColdHot',vmin=vmi,vmax=100)
axs[1].colorbar(m,label='%')
axs[1].format(boundinglat=50)
axs[4].format(boundinglat=-50)
axs[1].set_title('Fine fraction AOD')


m=axs[2].pcolormesh(lonsm,latsm,modis11,cmap='lajolla',vmin=0,vmax=0.15)
m1=axs[5].pcolormesh(lonsm,latsm,modis11,cmap='lajolla',vmin=0,vmax=0.15)
axs[2].format(boundinglat=50)
axs[5].format(boundinglat=-50)
axs[2].colorbar(m,label='AOD')
axs[2].set_title('Total AOD')

plt.show()
