########################
# Remy Lapere (03/11/2022)
# Figure 2/3 for Lapere et al., 2022 - CMIP6 SSaer
########################

import matplotlib.path as mpath
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,BoundaryNorm,LinearSegmentedColormap
import netCDF4
import numpy as np
import pandas as pd
import proplot as pplt
import cartopy.crs as ccrs
import cartopy
import matplotlib as mpl

mods = ['BCC-ESM','CESM','CNRM-ESM','EC-Earth','GISS','HadGEM','IPSL-CM6','MIROC-ES2L','MPI-ESM','MRI-ESM','NorESM','UKESM']

mdnm = ['mmrss_AERmon_BCC-ESM1_historical_r1i1p1f1_gn_185001-201412.nc','mmrss_AERmon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','mmrss_AERmon_CNRM-ESM2-1_historical_r1i1p1f2_gr_195001-201412.nc','mmrss_AERmon_EC-Earth3-AerChem_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_GISS-E2-1-H_historical_r1i1p3f1_gn_195101-201412.nc','mmrss_AERmon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_195001-201412.nc','mmrss_AERmon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gr_185001-201412.nc','mmrss_AERmon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','mmrss_AERmon_MPI-ESM-1-2-HAM_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_MRI-ESM2-0_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc']

# select the pole here
pole='Arctic'
#pole='Antarctic'

i=0
for md in mods:
    init = netCDF4.Dataset('../../../CMIP_DATA/mmrs/'+mdnm[i],'r')
    lons = init.variables['lon'][:]
    lats = init.variables['lat'][:]
    if md=='CESM':
        var = init.variables['mmrss'][-64*12:,-1,:,:]
    else:
        var = init.variables['mmrss'][-64*12:,0,:,:]
        
    outvar = np.mean(var,axis=0)*10**9
    
    if md=='IPSL-CM6':
        outvar = outvar/4.3 # IPSL-CM6 has wet SSaer mass -> need to convert

    if md=='CNRM-ESM':
        outvar=outvar/25.0 # to fit CNRM-ESM in the colorbar
        
    if pole=='Antarctic':
        vma=40
        ttickes = [0.01,0.05,0.1,0.5,1,5,10,13,16,19,22,25,28,31,34,37,40]
    else:
        vma=25
        ttickes = [0.01,0.05,0.1,0.5,1,5,8,10,13,16,19,22,25]
        
    outvar[np.isnan(outvar)] = 0.0
    
    colormap = plt.get_cmap('blues6')
    cols = [colormap(k) for k in np.linspace(0, 1, len(ttickes))]
    cols[3] = (0.6,0.6,0.6,1)
    cols[2] = (0.7,0.7,0.7,1)
    cols[1] = (0.8,0.8,0.8,1)
    cols[0] = (0.9,0.9,0.9,1)    
    cmap = LinearSegmentedColormap.from_list('testCmap', cols, N=len(ttickes))
    
    if pole=='Antarctic':
        ax1 = plt.subplot(3,4,i+1,projection=ccrs.SouthPolarStereo())
        ax1.set_extent([-180,180,-50,-90], ccrs.PlateCarree())
    else:
        ax1 = plt.subplot(3,4,i+1,projection=ccrs.NorthPolarStereo())
        ax1.set_extent([-180,180,55,90], ccrs.PlateCarree())    
        
    ax1.add_feature(cartopy.feature.LAND,edgecolor='k',zorder=4,facecolor=(0.5,0.5,0.5,0))
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax1.set_boundary(circle, transform=ax1.transAxes)

    m = ax1.pcolormesh(lons,lats,outvar,norm = BoundaryNorm(ttickes,len(ttickes)),transform=ccrs.PlateCarree(),cmap=cmap)
    ax1.set_title(md)
    
    i=i+1


plt.subplots_adjust(left=0.0,bottom=0.012,right=0.98,top=0.92,hspace=0.15,wspace=0)
plt.suptitle('Sea-salt aerosol surface mass mixing ratio - 1951-2014 average - '+pole)
plt.show()
