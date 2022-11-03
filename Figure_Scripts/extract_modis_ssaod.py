#!/bin/env python

import glob,os,sys
import matplotlib as mpl
mpl.use('tkagg')
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap,LogNorm
import pandas as pd
from matplotlib.collections import PatchCollection
from pyhdf.SD import SD,SDC
from matplotlib.patches import PathPatch
import scipy

def open_modis(fil):
    ff=SD(fil,SDC.READ)
    datasets_dic = ff.datasets()
    src=ff.select('AOD_550_Dark_Target_Deep_Blue_Combined_Mean_Mean')
    modis=src.get()
    for key, value in src.attributes().iteritems():
        if key == 'add_offset':
            add_offset = value  
        if key == 'scale_factor':
            scale_factor = value
    modis=modis*scale_factor+add_offset
    modis[modis<0.]=np.nan
    modis_lat=ff.select('YDim').get()
    modis_lon=ff.select('XDim').get()

    src_=ff.select('Aerosol_AE1_Ocean_JHisto_vs_Opt_Depth')
    modis_=src_.get()
    for key, value in src_.attributes().iteritems():
        if key == 'add_offset':
            add_offset_ = value
        if key == 'scale_factor':
            scale_factor_ = value
    modis_=modis_*scale_factor_+add_offset_
    modis_[modis_<0.]=np.nan
    
    src1=ff.select('AE1_Ocean_JHisto_Intervals')
    modis1=src1.get()
    modis_ = np.nansum(modis_[:,:,:,:],axis=0)
    modout = np.copy(modis_)*np.nan
    for h in range(len(modout[:,0,0])):
        modout[h,:,:] = modis_[h,:,:]*modis1[h]

    sp = np.nansum(modis_,axis=0)
    proba = modis_/sp
    modiss = np.nansum(modout,axis=0)/np.nansum(modis_,axis=0)

    return modis,modis_lat,modis_lon,modis1,proba,modiss


ddt = np.zeros(120)
star = 30

k=0
#for filepath in glob.iglob(r'../../../MYD08/CMIP/MYD08_M3.*20*.hdf'):
for filepath in glob.iglob(r'../../../MOD08/CMIP/MOD08_M3.*.20*.hdf'):
    ddt[k] = int(str(filepath)[star:star+7])# 6 char in the middle of the filename
    k=k+1

ddt = np.sort(ddt)
outmod = np.zeros((k,180,360))
ddts = np.zeros(k)
outmod_ = np.zeros((k,180,360))
prob = np.zeros((k,8,180,360))

j=0
for uu in range(120):
    #for filepath in glob.iglob(r'../../../MYD08/CMIP/MYD08_M3.A'+str(ddt[uu])+'*.hdf'):
    for filepath in glob.iglob(r'../../../MOD08/CMIP/MOD08_M3.A'+str(ddt[uu])+'*.hdf'):
        dd = str(filepath)[star:star+7]
        modis,lat,lon,aes,prob_,modiss = open_modis(filepath)
        outmod[j,:,:] = modis
        ddts[j] = int(dd)
        outmod_[j,:,:] = modiss
        prob[j,:,:,:] = prob_
        j=j+1

nc = netCDF4.Dataset('terra_ssaod.nc', 'w')
lat_dim = nc.createDimension('lat', len(lat))
lon_dim = nc.createDimension('lon', len(lon))
t_dim = nc.createDimension('time', np.shape(outmod)[0])
ae_dim = nc.createDimension('ae_bin', 8)
lat_var = nc.createVariable('lat', np.float64, ('lat'))
lat_var[:] = lat
lon_var = nc.createVariable('lon', np.float64, ('lon'))
lon_var[:] = lon
tnd = nc.createVariable('aodSS', np.float64, ('time','lat','lon'))
tnd[:,:,:] = outmod
times = nc.createVariable('time', np.float64, ('time'))
times[:] = ddts
aeprob = nc.createVariable('aeprob', np.float64, ('time','ae_bin','lat','lon'))
aeprob[:,:,:] = prob
aess = nc.createVariable('aes', np.float64, ('ae_bin'))
aess[:] = aes
ang = nc.createVariable('ang', np.float64, ('time','lat','lon'))
ang[:,:,:] = outmod_
nc.close()
