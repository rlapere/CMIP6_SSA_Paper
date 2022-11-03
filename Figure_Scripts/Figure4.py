########################
# Remy Lapere (03/11/2022)
# Figure 4 for Lapere et al., 2022 - CMIP6 SSaer
########################

from matplotlib.colors import LogNorm,BoundaryNorm,LinearSegmentedColormap
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy
import matplotlib as mpl
import matplotlib.style as style
import seaborn as sns
import proplot
from haversine import haversine, Unit
mpl.rcParams.update({'font.size': 12})

mods = ['BCC-ESM','CESM','CNRM-ESM','EC-Earth','GISS','HadGEM','IPSL-CM6','MIROC-ES2L','MPI-ESM','MRI-ESM','NorESM','UKESM']

mdnm = ['mmrss_AERmon_BCC-ESM1_historical_r1i1p1f1_gn_185001-201412.nc','mmrss_AERmon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','mmrss_AERmon_CNRM-ESM2-1_historical_r1i1p1f2_gr_195001-201412.nc','mmrss_AERmon_EC-Earth3-AerChem_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_GISS-E2-1-H_historical_r1i1p3f1_gn_195101-201412.nc','mmrss_AERmon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_195001-201412.nc','mmrss_AERmon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gr_185001-201412.nc','mmrss_AERmon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','mmrss_AERmon_MPI-ESM-1-2-HAM_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_MRI-ESM2-0_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc']

emiss = ['emiss_AERmon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','emiss_AERmon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','emiss_AERmon_CNRM-ESM2-1_historical_r1i1p1f2_gr_185001-201412.nc','emiss_AERmon_EC-Earth3-AerChem_historical_r1i1p1f1_gn_195001-201412.nc','emiss_AERmon_GISS-E2-1-H_historical_r1i1p3f1_gn_195101-201412.nc','emiss_AERmon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_195001-201412.nc','emiss_AERmon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gr_185001-201412.nc','emiss_AERmon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','emiss_AERmon_MPI-ESM-1-2-HAM_historical_r1i1p1f1_gn_195001-201412.nc','emiss_AERmon_MRI-ESM2-0_historical_r1i1p1f1_gn_185001-201412.nc','emiss_AERmon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc','emiss_AERmon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc']

aod = ['od550ss_AERmon_BCC-ESM1_historical_r1i1p1f1_gn_185001-201412.nc','od550ss_AERmon_CNRM-ESM2-1_historical_r1i1p1f2_gr_185001-201412.nc','od550ss_AERmon_EC-Earth3-AerChem_historical_r1i1p1f1_gn_195001-201412.nc','od550ss_AERmon_GISS-E2-1-H_historical_r1i1p3f1_gn_195101-201412.nc','od550ss_AERmon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gr_185001-201412.nc','od550ss_AERmon_MPI-ESM-1-2-HAM_historical_r3i1p1f1_gn_195001-201412.nc','od550ss_AERmon_MRI-ESM2-0_historical_r1i1p1f1_gn_185001-201412.nc','od550ss_AERmon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc']

modep = ['CESM','EC-Earth','GISS','IPSL-CM6','MIROC-ES2L','MPI-ESM','MRI-ESM','NorESM']

dryss = ['dryss_AERmon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','dryss_AERmon_EC-Earth3-AerChem_historical_r1i1p1f1_gn_195001-201412.nc','dryss_AERmon_GISS-E2-1-H_historical_r1i1p3f1_gn_195101-201412.nc','dryss_AERmon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gr_185001-201412.nc','dryss_AERmon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','dryss_AERmon_MPI-ESM-1-2-HAM_historical_r3i1p1f1_gn_195001-201412.nc','dryss_AERmon_MRI-ESM2-0_historical_r1i1p1f1_gn_185001-201412.nc','dryss_AERmon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc']

wetss = ['wetss_AERmon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','wetss_AERmon_EC-Earth3-AerChem_historical_r1i1p1f1_gn_195001-201412.nc','wetss_AERmon_GISS-E2-1-H_historical_r1i1p3f1_gn_195101-201412.nc','wetss_AERmon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gr_185001-201412.nc','wetss_AERmon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','wetss_AERmon_MPI-ESM-1-2-HAM_historical_r3i1p1f1_gn_195001-201412.nc','wetss_AERmon_MRI-ESM2-0_historical_r1i1p1f1_gn_185001-201412.nc','wetss_AERmon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc']

def gridtokm(lats,lons):
    surf = np.zeros((len(lats),len(lons)))
    lons_ = np.copy(lons)
    lons_[lons_>180] = lons_[lons_>180]-360.0
    for k in range(len(lats)-1):
        for l in range(len(lons_)-1):
            p1 = (lats[k],lons_[k])
            p2 = (lats[k],lons_[k+1])
            dx = haversine(p1,p2,unit=Unit.METERS)
            p1 = (lats[k],lons_[k])
            p2 = (lats[k+1],lons_[k])
            dy = haversine(p1,p2,unit=Unit.METERS)
            surf[k+1,l+1] = dx*dy
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

syms = ['^','o','2','*','<','>','x','+','d','s','v','1']

outm = np.zeros((len(mods),2))
oute=np.zeros((len(mods),2))
outod=np.zeros((len(mods),2))
outdep=np.zeros((len(mods),2))
outdepw=np.zeros((len(mods),2))

# pre-computed aerosol layer height
extn = pd.read_csv('extn.csv',header=None)[0].values
exts = pd.read_csv('exts.csv',header=None)[0].values

i=0
jj=0
kk=0
for md in mods:
    init = netCDF4.Dataset('../../../CMIP_DATA/mmrs/regridded_ocean/'+mdnm[i],'r')
    emi = netCDF4.Dataset('../../../CMIP_DATA/emiss/regridded_ocean/'+emiss[i],'r')
    if np.isin(md,['BCC-ESM','CNRM-ESM','EC-Earth','IPSL-CM6','MPI-ESM','MRI-ESM','NorESM','GISS']):
        od = netCDF4.Dataset('../../../CMIP_DATA/aod/regridded/'+aod[jj],'r')
        lons = od.variables['lon'][:]
        lats = od.variables['lat'][:]
        indn = np.where(lats>60)[0]
        inds = np.where(lats<-60)[0]
        ods = od.variables['od550ss'][-64*12:,:,:]
        outod[i,0] = np.mean(ods[:,indn,:])
        outod[i,1] = np.mean(ods[:,inds,:])
        jj=jj+1
    else:
        outod[i,0] = np.nan
        outod[i,1] = np.nan

    if np.isin(md,modep):
        dry = netCDF4.Dataset('../../../CMIP_DATA/dryss/'+dryss[kk],'r')
        wet = netCDF4.Dataset('../../../CMIP_DATA/wetss/'+wetss[kk],'r')
        llons_ = dry.variables['lon'][:]
        lats_ = dry.variables['lat'][:]
        indn_ = np.where(lats_>60)[0]
        inds_ = np.where(lats_<-60)[0]
        drys = dry.variables['dryss'][-64*12:,:,:]
        wets = wet.variables['wetss'][-64*12:,:,:]
        drys = np.mean(drys,axis=0)
        wets = np.mean(wets,axis=0)
        if np.min(drys)<0:
            drys = -drys
        if np.min(wets)<0:
            wets = -wets
        if md=='IPSL-CM6':
            drys = drys/4.3 # IPSL-CM6 provides SSaer at 80% RH instead for dry +> apply conversion factor
        gw_d = gridtokm(lats_,llons_)
        outdep[i,0] = np.nansum(drys[indn_,:]*gw_d[indn_,:])*3600*24*365/10**9
        outdep[i,1] = np.nansum(drys[inds_,:]*gw_d[inds_,:])*3600*24*365/10**9
        outdepw[i,0] = np.nansum(wets[indn_,:]*gw_d[indn_,:])*3600*24*365/10**9
        outdepw[i,1] = np.nansum(wets[inds_,:]*gw_d[inds_,:])*3600*24*365/10**9
        kk=kk+1
    else:
        outdep[i,0] = np.nan
        outdep[i,1] = np.nan
        outdepw[i,0] = np.nan
        outdepw[i,1] = np.nan

    lons = init.variables['lon'][:]
    lats = init.variables['lat'][:]
    lonse = emi.variables['lon'][:]
    latse = emi.variables['lat'][:]
    
    indn = np.where(lats>60)[0]
    inds = np.where(lats<-60)[0]
    indne = np.where(latse>60)[0]
    indse = np.where(latse<-60)[0]

    extn_ = extn[i]
    exts_ = exts[i]
    
    if md=='CESM':
        var = init.variables['mmrss'][-64*12:,-1,:,:]*10**9
    else:
        var = init.variables['mmrss'][-64*12:,0,:,:]*10**9

    resol = gridtokm(lats[indn],lons)
    resol1 = gridw(lats[indn],lons)
    resols = gridtokm(lats[inds],lons)
    resol1s = gridw(lats[inds],lons)
   
    emi_var = emi.variables['emiss'][-64*12:,:,:]
    
    if md=='IPSL-CM6':
        emi_var = emi_var/4.3
        var = var/4.3
        
    if md=='BCC-ESM':
         emi_var = emi_var*np.nan

    cmipemin = np.mean(emi_var[:,indn,:],axis=0)
    cmipemis = np.mean(emi_var[:,inds,:],axis=0)
    
    cmipn = np.nanmean(var[:,indn,:],axis=0)
    cmips = np.nanmean(var[:,inds,:],axis=0)

    cmipn[cmipemin==0.0]=np.nan
    cmips[cmipemis==0.0]=np.nan
    cmipn = np.nansum(cmipn*resol1)/np.sum(resol1)
    cmips = np.nansum(cmips*resol1s)/np.sum(resol1s)

    cmipn = np.nanmean(cmipn)*10
    cmips = np.nanmean(cmips)*10

    cmipemin = cmipemin*resol
    cmipemis = cmipemis*resols

    cmipemin = np.nansum(cmipemin)*3600.0*24.0*365.0/(10**9)
    cmipemis = np.nansum(cmipemis)*3600.0*24.0*365.0/(10**9)
    
    if md=='BCC-ESM':
        cmipemin = cmipemin*np.nan
        cmipemis = cmipemis*np.nan
            
    outm[i,0] = cmipemin
    outm[i,1] = cmipemis
    oute[i,0] = cmipn
    oute[i,1] = cmips
    i=i+1

# weighting for not taking CNRM-ESM into account in the colormap
ww = [1,1,0,1,1,1,1,1,1,1,1,1]

north = pd.DataFrame({'Mass emissions (Tg/yr)':outm[:,0]*ww/np.nanmax(outm[:,0]*ww),'Mass mixing ratio (x10 '+u"\u03bc"+'g/g)':oute[:,0]*ww/np.nanmax(oute[:,0]*ww),'Aerosol layer height (m)':extn/np.nanmax(extn),'AOD (x1000)':outod[:,0]/np.nanmax(outod[:,0]),'dep':outdep[:,0]/np.nanmax(outdep[:,0]),'depw':outdepw[:,0]/np.nanmax(outdepw[:,0])})
north1 = pd.DataFrame({'Mass emissions (Tg/yr)':outm[:,0].astype(int),'Mass mixing ratio ('+u"\u03bc"+'g/g)':oute[:,0].astype(int),'Aerosol layer height (m)':extn.astype(int),'AOD (x1000)':outod[:,0]*1000,'dep':outdep[:,0],'depw':outdepw[:,0]})

south = pd.DataFrame({'Mass emissions (Tg/yr)':outm[:,1]*ww/np.nanmax(outm[:,1]*ww),'Mass mixing ratio (x10 '+u"\u03bc"+'g/g)':oute[:,1]*ww/np.nanmax(oute[:,1]*ww),'Aerosol layer height (m)':exts/np.nanmax(exts),'AOD (x1000)':outod[:,1]/np.nanmax(outod[:,1]),'dep':outdep[:,1]/np.nanmax(outdep[:,1]),'depw':outdepw[:,1]/np.nanmax(outdepw[:,1])})
south1 = pd.DataFrame({'Mass emissions (Tg/yr)':outm[:,1].astype(int),'Mass mixing ratio ('+u"\u03bc"+'g/g)':oute[:,1].astype(int),'Aerosol layer height (m)':exts.astype(int),'AOD (x1000)':outod[:,1]*1000,'dep':outdep[:,1],'depw':outdepw[:,1]})

ttickes = np.linspace(0,1,100)
colormap = plt.get_cmap('blues6')
cols = [colormap(k) for k in np.linspace(0.0, 1, len(ttickes))]
labx = ['Mass emissions \n (Tg/yr)','Mass mixing \n ratio (x10 '+u"\u03bc"+'g/g)','Aerosol layer \n height (m)','AOD (x1000)','Dry deposition \n (Tg/year)','Wet deposition \n (Tg/year)']

cmap = LinearSegmentedColormap.from_list('testCmap', cols, N=len(ttickes))

ax=plt.subplot(131)
sns.heatmap(north,annot=north1,norm = BoundaryNorm(ttickes,len(ttickes)),cmap=cmap,fmt='.0f',linecolor='k',cbar=False,linewidths=0.5)
plt.yticks(np.arange(0.5,len(mods)+0.5,1),mods,rotation=0)
plt.xticks([0.5,1.5,2.5,3.5,4.5,5.5],labx,rotation=90)
ax.tick_params(which='minor',bottom=False,left=False)
ax.spines['bottom'].set_color('k')
ax.spines['top'].set_color('k')
ax.spines['left'].set_color('k')
ax.spines['right'].set_color('k')

ax.set_title('>60N')

ax2 = plt.subplot(132)
sns.heatmap(south,annot=south1,norm = BoundaryNorm(ttickes,len(ttickes)),cmap=cmap,fmt='.0f',linecolor='k',cbar=False,linewidths=0.5)
plt.yticks(np.arange(0.5,len(mods)+0.5,1),[],rotation=0)
plt.xticks([0.5,1.5,2.5,3.5,4.5,5.5],labx,rotation=90)
ax2.spines['bottom'].set_color('k')
ax2.spines['top'].set_color('k')
ax2.spines['left'].set_color('k')
ax2.spines['right'].set_color('k')
ax2.tick_params(left=False,which='both')
ax2.tick_params(which='minor',bottom=False)
ax2.set_title('<60S')

ax3 = plt.subplot(1,15,12)
lft = pd.read_csv('lftglobal.csv',header=None,names=['Lifetime (hour)'])
lft[lft['Lifetime (hour)']==0]=np.nan
ttickes = np.linspace(np.nanmin(lft['Lifetime (hour)'].values),np.nanmax(lft['Lifetime (hour)'].values),100)
sns.heatmap(lft,annot=lft,norm = BoundaryNorm(ttickes,len(ttickes)),cmap=cmap,fmt='.0f',linecolor='k',cbar=False,linewidths=0.5)
plt.yticks(np.arange(0.5,len(mods)+0.5,1),[],rotation=0)
plt.xticks([0.5],['Lifetime (hour)'],rotation=90)
ax3.spines['bottom'].set_color('k')
ax3.spines['top'].set_color('k')
ax3.spines['left'].set_color('k')
ax3.spines['right'].set_color('k')
ax3.set_title('Global')
ax3.tick_params(left=False,which='both')
ax3.tick_params(which='minor',bottom=False)
plt.show()
