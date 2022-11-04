###########
# Script to compute SSaer lifetime
# Lapere et al., 2022 - CMIP6 SSaer
##########

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,BoundaryNorm
import netCDF4
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy
import matplotlib as mpl
from haversine import haversine, Unit

mods = ['BCC-ESM','CESM','CNRM-ESM','EC-Earth','GISS','HadGEM','IPSL-CM6','MIROC-ES2L','MPI-ESM','MRI-ESM','NorESM','UKESM']

mdnm = ['mmrss_AERmon_BCC-ESM1_historical_r1i1p1f1_gn_185001-201412.nc','mmrss_AERmon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','mmrss_AERmon_CNRM-ESM2-1_historical_r1i1p1f2_gr_195001-201412.nc','mmrss_AERmon_EC-Earth3-AerChem_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_GISS-E2-1-H_historical_r1i1p3f1_gn_195101-201412.nc','mmrss_AERmon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_195001-201412.nc','mmrss_AERmon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gr_185001-201412.nc','mmrss_AERmon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','mmrss_AERmon_MPI-ESM-1-2-HAM_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_MRI-ESM2-0_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc']

emiss = ['emiss_AERmon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','emiss_AERmon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','emiss_AERmon_CNRM-ESM2-1_historical_r1i1p1f2_gr_185001-201412.nc','emiss_AERmon_EC-Earth3-AerChem_historical_r1i1p1f1_gn_195001-201412.nc','emiss_AERmon_GISS-E2-1-H_historical_r1i1p3f1_gn_195101-201412.nc','emiss_AERmon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_195001-201412.nc','emiss_AERmon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gr_185001-201412.nc','emiss_AERmon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','emiss_AERmon_MPI-ESM-1-2-HAM_historical_r1i1p1f1_gn_195001-201412.nc','emiss_AERmon_MRI-ESM2-0_historical_r1i1p1f1_gn_185001-201412.nc','emiss_AERmon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc','emiss_AERmon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc']

siconc = ['siconc_SImon_BCC-ESM1_historical_r1i1p1f1_gn_185001-201412.nc','siconc_SImon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','siconc_SImon_CNRM-ESM2-1_historical_r1i1p1f2_gn_185001-201412.nc','siconc_SImon_EC-Earth3-AerChem_historical_r1i1p1f1_gn_195001-201412.nc','siconc_SImon_GISS-E2-1-H_historical_r1i1p3f1_gr_195001-201412.nc','siconc_SImon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_195001-201412.nc','siconc_SImon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gn_185001-201412.nc','siconc_SImon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','siconc_SImon_MPI-ESM-1-2-HAM_historical_r1i1p1f1_gn_195001-201412.nc','siconc_SImon_MRI-ESM2-0_historical_r1i1p1f1_gn_185001-201412.nc','siconc_SImon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc','siconc_SImon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc']

def gridtokm(lats,lons):
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


def compute_lft(lats,lons):
    mm_ = mm[:,:,:]
    pm = ext_levs(md,init)
    pm = pm[:,:,:]
    
    if md=='CESM':
        dlev = pm[:-1,:,:]-pm[1:,:,:]
        load = mm_[1:,:,:]*dlev
        load = np.nansum(load,axis=0)
    else:
        dlev = pm[1:,:,:]-pm[:-1,:,:]
        load = mm_[:-1,:,:]*dlev
        load = np.nansum(load,axis=0)
    return load

lftn_ = np.zeros(len(mods))

i=0
for md in mods:
    init = netCDF4.Dataset('../../../CMIP_DATA/mmrs/regridded_ocean/'+mdnm[i],'r')
    emi = netCDF4.Dataset('../../../CMIP_DATA/emiss/regridded_ocean/'+emiss[i],'r')
    lons = init.variables['lon'][:]
    lats = init.variables['lat'][:]

    gridn = gridtokm(lats,lons)

    mm = init.variables['mmrss'][-12:,:,:,:]*1.2922
    mm = np.mean(mm,axis=0)
    if md=='IPSL-CM6':
        mm = mm/4.3
    
    levs = init.variables['lev'][:]
    emi_var = emi.variables['emiss'][-12:,:,:]
    emi_var = np.mean(emi_var,axis=0)
    
    loadn = compute_lft(lats,lons)
    loadn = loadn*gridn
    
    if md=='BCC-ESM':
         emi_var = emi_var*np.nan

    if md=='IPSL-CM6':
        emi_var = emi_var/4.3

    emi = emi_var*gridn
    
    if md!='BCC-ESM':
        lftn_[i] = np.nansum(loadn)/np.nansum(emi)/3600
 
    i = i+1

np.savetxt('global_SSaer_lifetime.csv',lftn_)

