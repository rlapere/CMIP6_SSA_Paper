########################
# Remy Lapere (03/11/2022)
# Figure 5 for Lapere et al., 2022 - CMIP6 SSaer
########################

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,BoundaryNorm
import netCDF4
import numpy as np
import proplot as pplt
import pandas as pd
import cartopy.crs as ccrs
import cartopy
import matplotlib as mpl
from haversine import haversine, Unit
from sklearn.metrics import explained_variance_score
mpl.rcParams.update({'font.size': 16})

mods = ['BCC-ESM','CESM','EC-Earth','GISS','HadGEM','IPSL-CM6','MIROC-ES2L','MPI-ESM','MRI-ESM','NorESM','UKESM']

mdnm = ['mmrss_AERmon_BCC-ESM1_historical_r1i1p1f1_gn_185001-201412.nc','mmrss_AERmon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','mmrss_AERmon_EC-Earth3-AerChem_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_GISS-E2-1-H_historical_r1i1p3f1_gn_195101-201412.nc','mmrss_AERmon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_195001-201412.nc','mmrss_AERmon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gr_185001-201412.nc','mmrss_AERmon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','mmrss_AERmon_MPI-ESM-1-2-HAM_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_MRI-ESM2-0_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc']

emiss = ['emiss_AERmon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','emiss_AERmon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','emiss_AERmon_EC-Earth3-AerChem_historical_r1i1p1f1_gn_195001-201412.nc','emiss_AERmon_GISS-E2-1-H_historical_r1i1p3f1_gn_195101-201412.nc','emiss_AERmon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_195001-201412.nc','emiss_AERmon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gr_185001-201412.nc','emiss_AERmon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','emiss_AERmon_MPI-ESM-1-2-HAM_historical_r1i1p1f1_gn_195001-201412.nc','emiss_AERmon_MRI-ESM2-0_historical_r1i1p1f1_gn_185001-201412.nc','emiss_AERmon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc','emiss_AERmon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc']

siconc = ['siconc_SImon_BCC-ESM1_historical_r1i1p1f1_gn_185001-201412.nc','siconc_SImon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','siconc_SImon_EC-Earth3-AerChem_historical_r1i1p1f1_gn_195001-201412.nc','siconc_SImon_GISS-E2-1-H_historical_r1i1p3f1_gr_195001-201412.nc','siconc_SImon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_195001-201412.nc','siconc_SImon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gn_185001-201412.nc','siconc_SImon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','siconc_SImon_MPI-ESM-1-2-HAM_historical_r1i1p1f1_gn_195001-201412.nc','siconc_SImon_MRI-ESM2-0_historical_r1i1p1f1_gn_185001-201412.nc','siconc_SImon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc','siconc_SImon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc']

sfcw=['sfcWind_Amon_BCC-ESM1_historical_r1i1p1f1_gn_185001-201412.nc','sfcWind_Amon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','sfcWind_Amon_EC-Earth3-AerChem_historical_r1i1p1f1_gr_195001-201412.nc','sfcWind_Amon_GISS-E2-1-H_historical_r1i1p3f1_gn_195101-201412.nc','sfcWind_Amon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_195001-201412.nc','sfcWind_Amon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gr_185001-201412.nc','sfcWind_Amon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','sfcWind_Amon_MPI-ESM-1-2-HAM_historical_r1i1p1f1_gn_195001-201412.nc','sfcWind_Amon_MRI-ESM2-0_historical_r1i1p1f1_gn_185001-201412.nc','sfcWind_Amon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc','sfcWind_Amon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc']

pblh = ['bldep_AERmon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','bldep_AERmon_EC-Earth3-AerChem_historical_r1i1p1f1_gn_195001-201412.nc','bldep_AERmon_GISS-E2-1-H_historical_r1i1p3f1_gn_195101-201412.nc','bldep_AERmon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_195001-201412.nc','bldep_AERmon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gr_185001-201412.nc','bldep_AERmon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','bldep_AERmon_MPI-ESM-1-2-HAM_historical_r1i1p1f1_gn_195001-201412.nc','bldep_AERmon_NorESM2-LM_historical_r2i1p1f1_gn_195001-201412.nc','bldep_AERmon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc']

def gridtokm(lats,lons):
    surf = np.zeros((12,len(lats),len(lons)))
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
            surf[:,k,l] = dx*dy
    surf[:,0] = surf[:,1]
    surf[0,:] = surf[1,:]
    return surf

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

f,axs = pplt.subplots(ncols=2,nrows=5,sharey=0)

axs0 = axs[0]
axs1 = axs[1]
axs2 = axs[2]
axs3 = axs[3]
axs4 = axs[4]
axs5 = axs[5]
axs6 = axs[6]
axs7 = axs[7]
axs8 = axs[8]
axs9 = axs[9]

syms = ['^','o','*','<','>','x','+','d','s','v','1']
colormap = plt.get_cmap('terrain')
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=[colormap(k) for k in np.linspace(0, 0.8, 10)])
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
col = [colors[0],colors[1],colors[2],colors[3],colors[4],colors[5],colors[6],colors[7],colors[8],colors[7],colors[9],colors[5]]

windout = np.zeros((11,12))
siout = np.zeros((11,12))
emiout = np.zeros((11,12))
windouts = np.zeros((11,12))
siouts = np.zeros((11,12))
emiouts = np.zeros((11,12))

i=0
kk=0
for md in mods:
    init = netCDF4.Dataset('../../../CMIP_DATA/mmrs/regridded_ocean/'+mdnm[i],'r')
    emi = netCDF4.Dataset('../../../CMIP_DATA/emiss/regridded_ocean/'+emiss[i],'r')
    si = netCDF4.Dataset('../../../CMIP_DATA/siconc/regridded_old/'+siconc[i],'r')
    ws = netCDF4.Dataset('../../../CMIP_DATA/sfcwind/regridded/'+sfcw[i],'r')

    lons = init.variables['lon'][:]
    lats = init.variables['lat'][:]
    lonse = emi.variables['lon'][:]
    latse = emi.variables['lat'][:]
    latsi = si.variables['lat'][:]
    lonsi = si.variables['lon'][:]
    latsw = ws.variables['lat'][:]
    lonsw = ws.variables['lon'][:]

    indn = np.where(lats>60)[0]
    inds = np.where(lats<-60)[0]
    indne = np.where(latse>60)[0]
    indse = np.where(latse<-60)[0]
    indni = np.where(latsi>60)[0]
    indsi = np.where(latsi<-60)[0]
    indnw = np.where(latsw>60)[0]
    indsw = np.where(latsw<-60)[0]


    gwn = gridw(lats[indn],lons)
    gws = gridw(lats[inds],lons)
    resol = gridtokm(lats,lons)
                    
    mm = init.variables['mmrss'][-64*12:,:,:,:]
    
    if md=='CESM':
        var = init.variables['mmrss'][-64*12:,-1,:,:]
    else:
        var = init.variables['mmrss'][-64*12:,0,:,:]

    emi_var = emi.variables['emiss'][-64*12:,:,:]

    var = var*10**9
    emi_var = emi_var

    emi_var = emi_var
    
    if md=='IPSL-CM6': # IPSL-CM6 provides SSaer at 80% RH instead for dry +> apply conversion factor
        emi_var = emi_var/4.3
        var = var/4.3
        
    si_var = si.variables['siconc'][-64*12:,:,:]
    si_var[si_var>100.0] = np.nan
    if md=='BCC-ESM':
        land = np.copy(si_var)
        land[~np.isnan(si_var)] = 1.0
    if md=='GISS':
        si_var = si_var*land
    si_var = 100.0-si_var

    sfcwind = ws.variables['sfcWind'][-64*12:,:,:]

    si_varn = si_var[:,indni,:]
    si_vars = si_var[:,indsi,:]
    
    windopocean = np.copy(sfcwind)
    windopocean[emi_var==0.0]=np.nan
    windopocean[np.isnan(emi_var)]=np.nan
    windopocean[np.isnan(si_var)]=np.nan    
    windopocean[si_var<10.0]=np.nan

    var[np.isnan(emi_var)] = np.nan
    var[emi_var==0.0] = np.nan
    var[np.isnan(si_var)]=np.nan
    var[si_var<10.0]=np.nan
    
    if np.isin(md,['CESM','EC-Earth','GISS','HadGEM','IPSL-CM6','MIROC-ES2L','MPI-ESM','NorESM','UKESM']):
        pb = netCDF4.Dataset('../../../CMIP_DATA/bldep/regridded_ocean/'+pblh[kk],'r')
        kk=kk+1
        pbl = pb.variables['bldep'][-64*12:,:,:]
        latp = pb.variables['lat'][:]
        indnp = np.where(latp>60)[0]
        indsp = np.where(latp<-60)[0]
        pbl_ = np.copy(pbl)
        pbl_[emi_var==0.0]=np.nan
        pbl_[np.isnan(emi_var)]=np.nan
        pbl_[np.isnan(si_var)]=np.nan
        pbl_[si_var<10.0]=np.nan
        cmippn = np.array([np.nanmean(pbl_[ii::12,indnp,:],axis=0) for ii in range(12)])
        cmipps = np.array([np.nanmean(pbl_[ii::12,indsp,:],axis=0) for ii in range(12)])
        cmippn = np.nanmean(cmippn*gwn,axis=(1,2))
        cmipps = np.nanmean(cmipps*gws,axis=(1,2))
    else:
        cmipps = np.zeros(12)*np.nan
        cmippn = cmipps
    
 
    cmipn = np.array([np.nanmean(var[i::12,indn,:],axis=0) for i in range(12)])
    cmipemin = np.array([np.nanmean(emi_var[i::12,indne,:],axis=0) for i in range(12)])
    cmipsin = np.array([np.nanmean(si_var[i::12,indni,:],axis=0) for i in range(12)])
    cmipwn = np.array([np.nanmean(windopocean[i::12,indnw,:],axis=0) for i in range(12)])

    cmips = np.array([np.nanmean(var[i::12,inds,:],axis=0) for i in range(12)])
    cmipemis = np.array([np.nanmean(emi_var[i::12,indse,:],axis=0) for i in range(12)])
    cmipsis = np.array([np.nanmean(si_var[i::12,indsi,:],axis=0) for i in range(12)])
    cmipws = np.array([np.nanmean(windopocean[i::12,indsw,:],axis=0) for i in range(12)])

    cmipn = np.nanmean(cmipn*gwn,axis=(1,2))
    cmips = np.nanmean(cmips*gws,axis=(1,2))
    cmipsin = np.nanmean(cmipsin*gwn,axis=(1,2))
    cmipsis = np.nanmean(cmipsis*gws,axis=(1,2))
    windn = np.nanmean(cmipwn*gwn,axis=(1,2))
    winds = np.nanmean(cmipws*gws,axis=(1,2))
    cmipemin = np.nansum(cmipemin*resol[:,indn,:],axis=(1,2))
    cmipemis = np.nansum(cmipemis*resol[:,inds,:],axis=(1,2))

    if md=='BCC-ESM':
        cmipemin = cmipemin*np.nan
        cmipemis = cmipemis*np.nan

    windout[i,:] = windn
    windouts[i,:] = winds
    emiout[i,:] = cmipemin
    emiouts[i,:] = cmipemis
    siout[i,:] = cmipsin
    siouts[i,:] = cmipsis

    cmap = plt.get_cmap("cividis")
    
    outout = cmipn
    outout = np.append((outout[0]+outout[-1])/2,outout)
    outout = np.append(outout,(outout[0]+outout[-1])/2)

    outouts = cmips
    outouts = np.append((outouts[0]+outouts[-1])/2,outouts)
    outouts = np.append(outouts,(outouts[0]+outouts[-1])/2)

    axs0.plot([0.5,1,2,3,4,5,6,7,8,9,10,11,12,12.5],outout,color=cmap((i+2)/13),marker=syms[i],label=md)
    axs1.plot([0.5,1,2,3,4,5,6,7,8,9,10,11,12,12.5],outouts,color=cmap((i+2)/13),marker=syms[i])
    
    outoute = cmipemin*3600*24/10**9 #Tg/month
    outoute = np.append((outoute[0]+outoute[-1])/2,outoute)
    outoute = np.append(outoute,(outoute[0]+outoute[-1])/2)

    outoutes = cmipemis*3600*24/10**9 #Tg/month
    outoutes = np.append((outoutes[0]+outoutes[-1])/2,outoutes)
    outoutes = np.append(outoutes,(outoutes[0]+outoutes[-1])/2)
        
    axs2.plot([0.5,1,2,3,4,5,6,7,8,9,10,11,12,12.5],outoute,color=cmap((i+2)/13),marker=syms[i])
    axs3.plot([0.5,1,2,3,4,5,6,7,8,9,10,11,12,12.5],outoutes,color=cmap((i+2)/13),marker=syms[i])

    outoute = cmipsin
    outoute = np.append((outoute[0]+outoute[-1])/2,outoute)
    outoute = np.append(outoute,(outoute[0]+outoute[-1])/2)

    outoutes =cmipsis
    outoutes = np.append((outoutes[0]+outoutes[-1])/2,outoutes)
    outoutes = np.append(outoutes,(outoutes[0]+outoutes[-1])/2)

    axs4.plot([0.5,1,2,3,4,5,6,7,8,9,10,11,12,12.5],outoute,color=cmap((i+2)/13),marker=syms[i])
    axs5.plot([0.5,1,2,3,4,5,6,7,8,9,10,11,12,12.5],outoutes,color=cmap((i+2)/13),marker=syms[i])

    outoute = windn
    outoute = np.append((outoute[0]+outoute[-1])/2,outoute)
    outoute = np.append(outoute,(outoute[0]+outoute[-1])/2)

    outoutes = winds
    outoutes = np.append((outoutes[0]+outoutes[-1])/2,outoutes)
    outoutes = np.append(outoutes,(outoutes[0]+outoutes[-1])/2)

    axs6.plot([0.5,1,2,3,4,5,6,7,8,9,10,11,12,12.5],outoute,color=cmap((i+2)/13),marker=syms[i])
    axs7.plot([0.5,1,2,3,4,5,6,7,8,9,10,11,12,12.5],outoutes,color=cmap((i+2)/13),marker=syms[i])

    outoutp = cmippn
    outoutp = np.append((outoutp[0]+outoutp[-1])/2,outoutp)
    outoutp = np.append(outoutp,(outoutp[0]+outoutp[-1])/2)

    outoutps = cmipps
    outoutps = np.append((outoutps[0]+outoutps[-1])/2,outoutps)
    outoutps = np.append(outoutps,(outoutps[0]+outoutps[-1])/2)

    axs8.plot([0.5,1,2,3,4,5,6,7,8,9,10,11,12,12.5],outoutp,color=cmap((i+2)/13),marker=syms[i])
    axs9.plot([0.5,1,2,3,4,5,6,7,8,9,10,11,12,12.5],outoutps,color=cmap((i+2)/13),marker=syms[i])
    
    i=i+1


siout = np.nanstd(siout,axis=0)
siouts = np.nanstd(siouts,axis=0)
windout = np.nanstd(windout,axis=0)
windouts = np.nanstd(windouts,axis=0)
emiout = np.nanstd(emiout,axis=0)
emiouts = np.nanstd(emiouts,axis=0)

axs0.grid(color='grey',axis='y')
axs1.grid(color='grey')
axs2.grid(color='grey')
axs3.grid(color='grey',axis='y')
axs4.grid(color='grey',axis='y')
axs5.grid(color='grey',axis='y')
axs6.grid(color='grey',axis='y')
axs7.grid(color='grey',axis='y')
axs8.grid(color='grey',axis='y')
axs9.grid(color='grey',axis='y')

axs0.set_xticks([])
axs1.set_xticks([])
axs2.set_xticks([])
axs3.set_xticks([])
axs4.set_xticks([])
axs5.set_xticks([])
axs6.set_xticks([])
axs7.set_xticks([])

axs0.set_xlim(0.5,12.5)
axs1.set_xlim(0.5,12.5)
axs2.set_xlim(0.5,12.5)
axs3.set_xlim(0.5,12.5)
axs4.set_xlim(0.5,12.5)
axs5.set_xlim(0.5,12.5)
axs6.set_xlim(0.5,12.5)
axs7.set_xlim(0.5,12.5)
axs8.set_xlim(0.5,12.5)
axs9.set_xlim(0.5,12.5)

axs0.set_title('a) Surface sea-salt mass mixing ratio $>$60N')
axs1.set_title('b) Surface sea-salt mass mixing ratio $<$60S')
axs2.set_title('c) Sea-salt mass emissions $>$60N')
axs3.set_title('d) Sea-salt mass emissions $<$60S')
axs4.set_title('e) Open ocean fraction $>$60N')
axs5.set_title('f) Open ocean fraction $<$60S')
axs6.set_title('g) Surface wind speed $>$60N')
axs7.set_title('h) Surface wind speed $<$60S')
axs8.set_title('i) Planetary boundary layer height $>$60N')
axs9.set_title('j) Planetary boundary layer height $<$60S')

axs8.set_xticks(np.arange(1,13))
axs9.set_xticks(np.arange(1,13))
axs8.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
axs9.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])

axs0.set_ylabel('Mass mixing ratio ('+u"\u03bc"+'g/g)')
axs2.set_ylabel('Emission rate (Tg/month)')
axs4.set_ylabel('Open ocean fraction (%)')
axs6.set_ylabel('Wind speed (m/s)')
axs8.set_ylabel('PBLH (m)')
f.legend(loc='bottom',ncol=6)
plt.suptitle('Arctic                                                                                  Antarctic')
plt.show()
