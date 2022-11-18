########################
# Remy Lapere (03/11/2022)
# Figure 11 and A7 for Lapere et al., 2022 - CMIP6 SSaer
########################

import pymannkendall as mk
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,BoundaryNorm,LinearSegmentedColormap
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
from haversine import haversine, Unit

mods = ['GISS',
        'HadGEM',
        'MIROC-ES2L',
        'MRI-ESM',
        'NorESM',
        'UKESM']

mdnmhist = ['mmrss_AERmon_GISS-E2-1-H_historical_r1i1p3f1_gn_195101-201412.nc',
            'mmrss_AERmon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_195001-201412.nc',
            'mmrss_AERmon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc',
            'mmrss_AERmon_MRI-ESM2-0_historical_r1i1p1f1_gn_195001-201412.nc',
            'mmrss_AERmon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc',
            'mmrss_AERmon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc']

mdnm585 = ['mmrss_AERmon_GISS-E2-1-H_ssp585_r1i1p3f1_gn_201501-210012.nc',
           'mmrss_AERmon_HadGEM3-GC31-LL_ssp585_r1i1p1f3_gn_201501-210012.nc',
           'mmrss_AERmon_MIROC-ES2L_ssp585_r1i1p1f2_gn_201501-210012.nc',
           'mmrss_AERmon_MRI-ESM2-0_ssp585_r1i1p1f1_gn_201501-210012.nc',
           'mmrss_AERmon_NorESM2-LM_ssp585_r1i1p1f1_gn_201501-210012.nc',
           'mmrss_AERmon_UKESM1-0-LL_ssp585_r1i1p1f2_gn_201501-210012.nc']

mdnm126 = ['mmrss_AERmon_GISS-E2-1-H_ssp126_r1i1p3f1_gn_201501-210012.nc',
           'mmrss_AERmon_HadGEM3-GC31-LL_ssp126_r1i1p1f3_gn_201501-210012.nc',
           'mmrss_AERmon_MIROC-ES2L_ssp126_r1i1p1f2_gn_201501-210012.nc',
           'mmrss_AERmon_MRI-ESM2-0_ssp126_r1i1p1f1_gn_201501-210012.nc',
           'mmrss_AERmon_NorESM2-LM_ssp126_r1i1p1f1_gn_201501-210012.nc',
           'mmrss_AERmon_UKESM1-0-LL_ssp126_r1i1p1f2_gn_201501-210012.nc']

def filt(da):
    return savgol_filter(da,19,3)
    

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


mmrs_ = netCDF4.Dataset('../../../CMIP_DATA/mmrs/regridded/'+mdnmhist[0],'r')
lons = mmrs_.variables['lon'][:]
lats = mmrs_.variables['lat'][:]

si0_ = netCDF4.Dataset('trend_siconc.nc','r')
lonssi = si0_.variables['lon'][:]
latssi = si0_.variables['lat'][:]
trendssi = si0_.variables['trend'][:,:]

emi0_ = netCDF4.Dataset('trend_emiss_rel_six_indiv.nc','r')
lonemi = emi0_.variables['lon'][:]
latemi = emi0_.variables['lat'][:]
trendemi = emi0_.variables['trend'][:,:]*12*100*10

w0_ = netCDF4.Dataset('trend_sfcwind.nc','r')
lonw = w0_.variables['lon'][:]
latw = w0_.variables['lat'][:]
trendw = w0_.variables['trend'][:,:]


fig, axs = plot.subplots(ncols=2,nrows=1,proj={1:'npaeqd',2:'spaeqd'},sharey=False,sharex=False)
axs.format(suptitle='Historical trend in sea-salt aerosol mass emission flux',coast=True, latlines=10,coastcolor='black',land=True,landcolor='white')

indn = np.where(latemi>50)[0]
inds = np.where(latemi<-60)[0]

vman = np.nanpercentile(np.abs(trendemi[indn,:].flatten()),99)
vmas = np.nanpercentile(np.abs(trendemi[inds,:].flatten()),99)

a=axs[0].pcolormesh(lonemi,latemi,trendemi,cmap='ColdHot',vmin=-vmas,vmax=vmas,extend='both')
axs[0].format(boundinglat=60,title='Arctic 1951$-$2014')
axs[0].colorbar(a,shrink=0.8,label='%$\,$decade$^{-1}$')

b=axs[1].pcolormesh(lonemi,latemi,trendemi,cmap='ColdHot',vmin=-vmas,vmax=vmas,extend='both')
axs[1].format(boundinglat=-60,title='Antarctic 1951$-$2014')

plt.show()


mmr585 = netCDF4.Dataset('trend_mmrss_rel_ssp585_indiv.nc','r')
lonm = mmr585.variables['lon'][:]
latm = mmr585.variables['lat'][:]
trendm = mmr585.variables['trend'][:,:]*12*100*10

fig, axs = plot.subplots(ncols=1,nrows=2,proj={1:'robin',2:'robin'})
axs.format(suptitle='Future trend in sea salt aerosol and sea-ice concentration',coast=True, latlines=10,coastcolor='black')

vman = np.nanpercentile(np.abs(trendm.flatten()),99)

a=axs[0].pcolormesh(lonm,latm,trendm,cmap='ColdHot',vmin=-vman,vmax=vman,extend='both',N=20)
axs[0].format(title='')
axs[0].colorbar(a,shrink=0.8,label='%$\,$decade$^{-1}$')

si0_ = netCDF4.Dataset('trend_mmrss_rel_ssp585_siconc.nc','r')
lonssi = si0_.variables['lon'][:]
latssi = si0_.variables['lat'][:]
trendssi = si0_.variables['trend'][:,:]*12*10
vmasi = np.nanmax(np.abs(trendssi.flatten()))

colormap = plt.get_cmap('blues6')
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=[colormap(k) for k in np.linspace(0, 1, 10)])
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
cmap = LinearSegmentedColormap.from_list('testCmap', np.flip(colors), N=10)

a=axs[1].pcolormesh(lonssi,latssi,trendssi,cmap='blues6_r',vmin=-vmasi,vmax=0)
axs[1].colorbar(a,shrink=0.8,label='Sea-ice concentration - %$\,$decade$^{-1}$')

plt.show()

fig, axs = plot.subplots(ncols=2,nrows=1,proj={1:'cartesian',2:'cartesian'},sharey=False,sharex=False)
axs.format(suptitle='Historical and future trend in sea salt aerosol mass mixing ratio',coast=True,coastcolor='black',land=True,landcolor='white')

ensn26 = np.zeros((len(mods),150))
enss26 = np.zeros((len(mods),150))
ensn85 = np.zeros((len(mods),150))
enss85 = np.zeros((len(mods),150))
ensms85 = np.zeros((len(mods),150))
ensmn85 = np.zeros((len(mods),150))

i=0
for md in mods:
    s85 = netCDF4.Dataset('../../../CMIP_DATA/ssp/'+mdnm585[i],'r')
    s26 = netCDF4.Dataset('../../../CMIP_DATA/ssp/'+mdnm126[i],'r')

    lons85 = s85.variables['lon'][:]
    lats85 = s85.variables['lat'][:]
    
    hist = netCDF4.Dataset('../../../CMIP_DATA/mmrs/'+mdnmhist[i],'r')
    lons = hist.variables['lon'][:]
    lats = hist.variables['lat'][:]

    varhist = hist.variables['mmrss'][-64*12:,0,:,:]
    var85 = s85.variables['mmrss'][:,0,:,:]
    var26 = s26.variables['mmrss'][:,0,:,:]
    
    indn = np.where(lats>60)[0]
    inds = np.where(lats<-60)[0]
    indn85 = np.where(lats85>60)[0]
    inds85 = np.where(lats85<-60)[0]

    gridws = gridw(lats[inds],lons)
    gridwn = gridw(lats[indn],lons)

    n85 = np.mean(var85[:,indn85,:]*gridwn,axis=(1,2))
    s85 = np.mean(var85[:,inds85,:]*gridws,axis=(1,2))

    n26 = np.mean(var26[:,indn85,:]*gridwn,axis=(1,2))
    s26 = np.mean(var26[:,inds85,:]*gridws,axis=(1,2))

    nhist = np.mean(varhist[:,indn,:]*gridwn,axis=(1,2))
    shist = np.mean(varhist[:,inds,:]*gridws,axis=(1,2))
        
    df = pd.DataFrame()
    df['date'] = pd.date_range('1951-01-01','2100-12-31',freq='m')
    df['year'] = df.date.dt.year
    df['month'] = df.date.dt.month
    df['n85'] = np.append(nhist,n85)
    df['s85'] = np.append(shist,s85)
    df['n26'] = np.append(nhist,n26)
    df['s26'] = np.append(shist,s26)

    cmipn85 = df.groupby('year',as_index=False).mean()['n85'].values
    cmips85 = df.groupby('year',as_index=False).mean()['s85'].values
    cmipn26 = df.groupby('year',as_index=False).mean()['n26'].values
    cmips26 = df.groupby('year',as_index=False).mean()['s26'].values
    
    cmipn85 = cmipn85/np.mean(cmipn85[:20])
    cmips85 = cmips85/np.mean(cmips85[:20])

    cmipn26 = cmipn26/np.mean(cmipn26[:20])
    cmips26 = cmips26/np.mean(cmips26[:20])

    cmipn26 = filt(cmipn26)
    cmips26 = filt(cmips26)
    cmipn85 = filt(cmipn85)
    cmips85 = filt(cmips85)

    ensn26[i,:] = cmipn26
    enss26[i,:] = cmips26
    ensn85[i,:] = cmipn85
    enss85[i,:] = cmips85

    lb = '_no_legend_'
    lb1 = '_no_legend_'
    lb0 = '_no_legend_'
        
    axs[0].plot(np.arange(63,150),cmipn85[63:],label=lb,color='pink',linewidth=1,alpha=0.6)
    axs[1].plot(np.arange(63,150),cmips85[63:],color='pink',linewidth=1,alpha=0.6)

    axs[0].plot(range(64),cmipn85[:64],label=lb0,color='gray',linewidth=1,alpha=1)
    axs[1].plot(range(64),cmips85[:64],color='gray',linewidth=1,alpha=1)
    
    axs[0].plot(np.arange(63,150),cmipn26[63:],label=lb1,color='lightblue',linewidth=1,alpha=0.6)
    axs[1].plot(np.arange(63,150),cmips26[63:],color='lightblue',linewidth=1,alpha=0.6)

    i=i+1

axs[0].plot(np.arange(63,150),np.mean(ensn85[:,63:],axis=0),label='SSP585',color='cherry',linewidth=3,alpha=1)
axs[1].plot(np.arange(63,150),np.mean(enss85[:,63:],axis=0),color='cherry',linewidth=3,alpha=1)

axs[0].plot(range(64),np.mean(ensn85[:,:64],axis=0),label='Historical',color='black',linewidth=3,alpha=1)
axs[1].plot(range(64),np.mean(enss85[:,:64],axis=0),color='black',linewidth=3,alpha=1)

axs[0].plot(np.arange(63,150),np.mean(ensn26[:,63:],axis=0),label='SSP126',color='ocean',linewidth=3,alpha=1)
axs[1].plot(np.arange(63,150),np.mean(enss26[:,63:],axis=0),color='ocean',linewidth=3,alpha=1)

axs[0].legend(ncol=1)
axs[0].grid(b=True,axis='both',color='gray')
axs[0].set_xticks(range(150)[::20])
axs[0].set_xticklabels(np.arange(1951,2101,20).astype(str))
axs[0].format(title='North of 60$^{\circ}$N average',ylabel='Relative evolution wrt 1951$-$1971')
axs[0].set_xlabel('Year')

axs[1].set_xticks(range(150)[::20])
axs[1].set_xticklabels(np.arange(1951,2101,20).astype(str))
axs[1].grid(b=True,axis='both',color='gray')
axs[1].format(title='South of 60$^{\circ}$S average')
axs[1].set_xlabel('Year')

plt.show()
