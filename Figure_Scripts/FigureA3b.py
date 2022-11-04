########################
# Remy Lapere (03/11/2022)
# Figure A3 right panel for Lapere et al., 2022 - CMIP6 SSaer
########################

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,BoundaryNorm
import netCDF4
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy
import matplotlib as mpl
mpl.rcParams.update({'font.size': 16})
from haversine import haversine, Unit

def gridw(lats,lons):
    surf = np.zeros((len(lats),len(lons)))
    lons_ = np.copy(lons)
    lons_[lons_>180] = lons_[lons_>180]-360.0
    for k in range(len(lats)-1):
        for l in range(len(lons)-1):
            p1 = (lats[k],lons_[k])
            p2 = (lats[k],lons_[k+1])
            dx = haversine(p1,p2,unit=Unit.METERS)
            p1 = (lats[k],lons_[k])
            p2 = (lats[k+1],lons_[k])
            dy = haversine(p1,p2,unit=Unit.METERS)
            surf[k,l] = dx*dy
    surf = surf/np.max(surf)
    if np.max(lats)<0:
        surf[-1,:]=1
    if np.min(lats)>0:
        surf[-1,:]=np.min(surf[surf>0])
    return surf


def ext_levs(md,init):
    if np.isin(md,['MPI-ESM','EC-Earth','BCC-ESM','GISS','MIROC-ES2L','MRI-ESM','NorESM','IPSL-CM6','CNRM-ESM']):
        levs = init.variables['lev'][:]
        if md=='IPSL-CM6':
            init0 = netCDF4.Dataset('../../../CMIP_DATA/mmrs/mmrss_AERmon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gr_185001-201412.nc','r')
            a = init0.variables['ap'][:]
            b = init0.variables['b'][:]
        else:
            a = init.variables['ap'][:]
            b = init.variables['b'][:]
        ps = init.variables['ps'][-65*12:,:,:]
        ps = np.mean(ps,axis=0)
        b_tile = np.tile(b,(np.shape(ps)[0],np.shape(ps)[1],1))
        ps_tile = np.repeat(ps[:,:,np.newaxis],len(b),axis=2)
        p = a+b_tile*ps_tile
        pm = np.zeros((levs.shape[0],ps.shape[0],ps.shape[1]))
        for ll in range(len(levs)):
            pm[ll,:,:] = 44330*(1-(p[:,:,ll]/101325)**(1/5.255))
    elif md=='CESM':
        levloc = netCDF4.Dataset('../../../CMIP_DATA/mmrs/mmrss_AERmon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','r')
        levs = init.variables['lev'][:]
        ps = init.variables['ps'][-65*12:,:,:]
        a = levloc.variables['a'][:]
        b = levloc.variables['b'][:]
        ps = np.mean(ps,axis=0)
        b_tile = np.tile(b,(np.shape(ps)[0],np.shape(ps)[1],1))
        ps_tile = np.repeat(ps[:,:,np.newaxis],len(b),axis=2)
        p = a+b_tile*ps_tile
        pm = np.zeros((levs.shape[0],ps.shape[0],ps.shape[1]))
        for ll in range(len(levs)):
            pm[ll,:,:] = 44330*(1-(p[:,:,ll]/101325)**(1/5.255))
    elif np.isin(md,['UKESM','HadGEM']):
        levs = init.variables['lev'][:]
        pm = np.zeros((levs.shape[0],init.variables['lat'].shape[0],init.variables['lon'].shape[0]))
        for ll in range(len(levs)):
            pm[ll,:,:] = levs[ll]
    return pm




mods = ['BCC-ESM','CESM','EC-Earth','GISS','HadGEM','IPSL-CM6','MIROC-ES2L','MPI-ESM','MRI-ESM','NorESM','UKESM']

mdnm = ['cl_Amon_BCC-ESM1_historical_r1i1p1f1_gn_185001-201412.nc','cl_Amon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','cl_Amon_EC-Earth3-AerChem_historical_r1i1p1f1_gr_195001-201412.nc','cl_Amon_GISS-E2-1-H_historical_r1i1p1f1_gn_195101-201412.nc','cl_Amon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_195001-201412.nc','cl_Amon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gr_185001-201412.nc','cl_Amon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','cl_Amon_MPI-ESM-1-2-HAM_historical_r1i1p1f1_gn_195001-201412.nc','cl_Amon_MRI-ESM2-0_historical_r1i1p1f1_gn_195001-201412.nc','cl_Amon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc','cl_Amon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc']

sic = ['siconc_SImon_BCC-ESM1_historical_r1i1p1f1_gn_185001-201412.nc','siconc_SImon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','siconc_SImon_EC-Earth3-AerChem_historical_r1i1p1f1_gn_195001-201412.nc','siconc_SImon_GISS-E2-1-H_historical_r1i1p3f1_gr_195001-201412.nc','siconc_SImon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_195001-201412.nc','siconc_SImon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gn_185001-201412.nc','siconc_SImon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','siconc_SImon_MPI-ESM-1-2-HAM_historical_r1i1p1f1_gn_195001-201412.nc','siconc_SImon_MRI-ESM2-0_historical_r1i1p1f1_gn_185001-201412.nc','siconc_SImon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc','siconc_SImon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc']

tarlev = np.logspace(1,4,100)

mrnhivout = np.zeros((len(tarlev),2))
mrnhivout[:,0] = 1000
mrshivout = np.zeros((len(tarlev),2))
mrshivout[:,0] = 1000
mrneteout = np.zeros((len(tarlev),2))
mrneteout[:,0] = 1000
mrseteout = np.zeros((len(tarlev),2))
mrseteout[:,0] = 1000

ax=plt.subplot(121)
ax1=plt.subplot(122)

i=0
for md in mods:
    pre = netCDF4.Dataset('../../../CMIP_DATA/mmrs/regridded_ocean/'+mdnm[i],'r')
    si = netCDF4.Dataset('../../../CMIP_DATA/siconc/regridded_old/siconc_SImon_BCC-ESM1_historical_r1i1p1f1_gn_185001-201412.nc','r')
    siconc = np.nanmean(si.variables['siconc'][:],axis=0)
    masko = np.copy(siconc)
    masko[masko>100] = np.nan
    masko[masko>50] = np.nan
    masko[~np.isnan(masko)] = 1.0
    lons = pre.variables['lon'][:]
    lats = pre.variables['lat'][:]
    indn = np.where(lats>60)[0]
    inds = np.where(lats<-60)[0]
    gwn = gridw(lats[indn],lons)
    gws = gridw(lats[inds],lons)
    mpr = pre.variables['mmrss'][-65*12:,:,:,:]
    if md=='IPSL-CM6':
        mpr = mpr/4.3
    for hh in range(np.shape(mpr)[1]):
        mpr[:,hh,:,:] = mpr[:,hh,:,:]*masko
    seas = range(len(mpr[:,0,0,0]))
    hiv = np.where((np.mod(seas,12)==0) | (np.mod(seas,12)==1) | (np.mod(seas,12)==11))[0]
    ete = np.where((np.mod(seas,12)==5) | (np.mod(seas,12)==6) | (np.mod(seas,12)==7))[0]
    spr = np.where((np.mod(seas,12)==2) | (np.mod(seas,12)==3) | (np.mod(seas,12)==4))[0]
    fal = np.where((np.mod(seas,12)==8) | (np.mod(seas,12)==9) | (np.mod(seas,12)==10))[0]
    sea = hiv
    mprhiv = mpr[hiv,:,:,:]
    mprete = mpr[ete,:,:,:]
    mrsete=np.nanmean(mprete[:,:,inds,:]*gws,axis=(0,2,3))
    mrnete=np.nanmean(mprete[:,:,indn,:]*gwn,axis=(0,2,3))
    mrshiv=np.nanmean(mprhiv[:,:,inds,:]*gws,axis=(0,2,3))
    mrnhiv=np.nanmean(mprhiv[:,:,indn,:]*gwn,axis=(0,2,3))
    
    levs = ext_levs(md,pre)
    for hh in range(np.shape(mpr)[1]):
        levs[hh,:,:] = levs[hh,:,:]*masko
    ml = np.nanmean(levs[:,indn,:]*gwn,axis=(1,2))
    mls = np.nanmean(levs[:,inds,:]*gws,axis=(1,2))
    
    if md=='CESM':
        ml = np.flip(ml)
        mrnhiv = np.flip(mrnhiv)
        mls = np.flip(mls)
        mrshiv = np.flip(mrshiv)
        mrnete = np.flip(mrnete)
        mrsete = np.flip(mrsete)
        
    mrnhivint = np.interp(tarlev,ml-np.min(ml),mrnhiv)
    mrshivint = np.interp(tarlev,mls-np.min(mls),mrshiv)
    mrneteint = np.interp(tarlev,ml-np.min(ml),mrnete)
    mrseteint = np.interp(tarlev,mls-np.min(mls),mrsete)

    mrnhivout[:,0] = np.minimum(mrnhivout[:,0],mrnhivint)
    mrnhivout[:,1] = np.maximum(mrnhivout[:,1],mrnhivint)
    mrshivout[:,0] = np.minimum(mrshivout[:,0],mrshivint)
    mrshivout[:,1] = np.maximum(mrshivout[:,1],mrshivint)
    mrneteout[:,0] = np.minimum(mrneteout[:,0],mrneteint)
    mrneteout[:,1] = np.maximum(mrneteout[:,1],mrneteint)
    mrseteout[:,0] = np.minimum(mrseteout[:,0],mrseteint)
    mrseteout[:,1] = np.maximum(mrseteout[:,1],mrseteint)
  

    ax.plot(mrnhivint,tarlev,color='teal',alpha=1,zorder=1)
    ax.plot(mrneteint,tarlev,color='darkgoldenrod',alpha=1,zorder=2)
    ax.set_ylim(10,10000)
    ax.set_yscale('log')
    ax.set_xlim(0,60)

    if i==0:
        lab = 'Jun'+'\u2013'+'Aug'
        lab1= 'Dec'+'\u2013'+'Feb'
    else:
        lab='_no_legend_'
        lab1 = '_no_legend_'
    
    ax1.plot(mrseteint,tarlev,label=lab,color='darkgoldenrod',alpha=1,zorder=2)
    ax1.plot(mrshivint,tarlev,label=lab1,color='teal',alpha=1,zorder=1)
    ax1.set_ylim(10,10000)
    ax1.set_yscale('log')
    ax1.set_xlim(0,60)
    
    i = i+1


ax.fill_betweenx(tarlev,mrnhivout[:,0],mrnhivout[:,1],alpha=0.25,color='teal',zorder=0)
ax.fill_betweenx(tarlev,mrneteout[:,0],mrneteout[:,1],alpha=0.5,color='darkgoldenrod',zorder=1)

ax1.fill_betweenx(tarlev,mrshivout[:,0],mrshivout[:,1],alpha=0.25,color='teal',zorder=0)
ax1.fill_betweenx(tarlev,mrseteout[:,0],mrseteout[:,1],alpha=0.5,color='darkgoldenrod',zorder=1)

ax1.set_title('<60S')
ax.set_title('>60N')

ax.set_ylabel('Altitude (m)')
ax.set_xlabel('Cloud fraction (%)')
ax1.set_xlabel('Cloud fraction (%)')

plt.legend()
plt.show()
