##########
# Remy Lapere
# Example script for computing trends in Lapere et al., 2022 - CMIP6 SSaer
##########

import pymannkendall as mk
import netCDF4
import numpy as np
import pandas as pd

mods = ['GISS','HadGEM','MIROC-ES2L','MRI-ESM','NorESM','UKESM']

mdnmhist = ['mmrss_AERmon_GISS-E2-1-H_historical_r1i1p3f1_gn_195101-201412.nc',
            'mmrss_AERmon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_195001-201412.nc',
            'mmrss_AERmon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc',
            'mmrss_AERmon_MRI-ESM2-0_historical_r1i1p1f1_gn_195001-201412.nc',
            'mmrss_AERmon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc',
            'mmrss_AERmon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc']

def comp_trend(din,ll,lo):
    pmk = mk.original_test(din[:,ll,lo])
    if pmk.h==False:
        return np.nan
    else:
        return pmk.slope

mmrs_ = netCDF4.Dataset('../../../CMIP_DATA/mmrs/regridded/'+mdnmhist[0],'r')
lons = mmrs_.variables['lon'][:]
lats = mmrs_.variables['lat'][:]
mmrs = mmrs_.variables['mmrss'][-64*12:,:]
vv_mm = np.zeros((len(mods),np.shape(mmrs)[0],len(lats),len(lons)))

# store the ensemble
i=0
for md in mods:
    mmrs_ = netCDF4.Dataset('../../../CMIP_DATA/mmrs/regridded/'+mdnmhist[i],'r')
    mmrs = mmrs_.variables['mmrss'][:,0,:,:]
    vv_mm[i,:,:,:] = mmrs
    i=i+1

ensemble = np.nanmean(vv_mm,axis=0)

trends = np.zeros((len(lats),len(lons)))
for la in range(len(lats)):
    for lo in range(len(lons)):
        trends[la,lo] = comp_trend(ensemble,la,lo)

nc = netCDF4.Dataset('trend.nc', 'w')
lat_dim = nc.createDimension('lat', len(lats))
lon_dim = nc.createDimension('lon', len(lons))
lat_var = nc.createVariable('lat', np.float64, ('lat'))
lat_var[:] = lats
lon_var = nc.createVariable('lon', np.float64, ('lon'))
lon_var[:] = lons
tnd = nc.createVariable('trend', np.float64, ('lat','lon'))
tnd[:,:] = trends
nc.close()
