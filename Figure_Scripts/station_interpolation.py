########################
# Remy Lapere (03/11/2022)
# extraction of SSAer station annual cycles in CMIP6 models for Lapere et al., 2022 - CMIP6 SSaer
########################

import numpy as np
import pandas as pd
import netCDF4

var = 'mmrss'

lgd = ['BCC','CESM','EC','GISS','HADGEM','IPSL','MIROC','MPI','MRI','NORESM','UKESM1']
lgd1 = ['BCC-ESM','CESM','EC-Earth','GISS','HadGEM','IPSL-CM6','MIROC-ES2L','MPI-ESM','MRI-ESM','NorESM','UKESM','ENS','CAMS','MERRA2']


models = ['mmrss_AERmon_BCC-ESM1_historical_r1i1p1f1_gn_185001-201412.nc','mmrss_AERmon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','mmrss_AERmon_EC-Earth3-AerChem_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_GISS-E2-1-H_historical_r1i1p3f1_gn_195101-201412.nc','mmrss_AERmon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_195001-201412.nc','mmrss_AERmon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gr_185001-201412.nc','mmrss_AERmon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','mmrss_AERmon_MPI-ESM-1-2-HAM_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_MRI-ESM2-0_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc']

ind=1

def find_point(la,lo,lats_,lons_):
    la_ = np.argmin(np.abs(lats_-la))
    lo_ = np.argmin(np.abs(lons_-lo))
    return la_,lo_

for st in stats:
    i=0
    for md in models:
        init = netCDF4.Dataset('../../../CMIP_DATA/mmrs/'+md,'r')
        laloc = init.variables['lat'][:]
        loloc = init.variables['lon'][:]
        ssa = init.variables['mmrss'][:,0,:,:]*10**9
        if lgd[i]=='CESM':
            ssa=init.variables['mmrss'][:,-1,:,:]*10**9
        if (st[2]<0.0) and np.max(loloc)>190.0:
            stloc = st[2]+360
        else:
            stloc = st[2]
        x,y = find_point(st[1],stloc,laloc,loloc)
        ssa = ssa[-65*12:,x,y]
        cmip = np.array([np.mean(ssa[mm::12]) for mm in range(12)])
        df_cmip = pd.DataFrame({'ssa':cmip})
        df_cmip.to_csv(st[0]+'_'+lgd1[i]+'.csv')
        i=i+1
