import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import proplot as plot
import pandas as pd
from matplotlib.colors import LogNorm
from scipy.stats import pearsonr
import matplotlib as mpl

def open_modis(cas):
    if cas=='aqua':
        fld = 'aqua_ssaod.nc'
        modreg = 'aqua_ssaod_reg.nc'#CMIP6 re-gridded MODIS
    elif cas=='terra':
        fld = 'terra_ssaod.nc'
        modreg = 'terra_ssaod_reg.nc'
 
    md = netCDF4.Dataset(fld,'r')
    modis = md.variables['aodSS'][:]
    latsm = md.variables['lat'][:]
    lonsm = md.variables['lon'][:]
    ddts = md.variables['time'][:]
    aeg = md.variables['ang'][:]
    poids = md.variables['aeprob'][:]
    aes = md.variables['aes'][:]    
    
    mdreg = netCDF4.Dataset(modreg,'r')
    modisr = mdreg.variables['aodSS'][:]
    latsmr = mdreg.variables['lat'][:]
    lonsmr = mdreg.variables['lon'][:]

    modmask = np.copy(modisr)
    modmask[:,np.abs(latsmr)>80,:] = np.nan
    modmask[~np.isnan(modmask)] = 1.0

    thresh = 2
    poids1 = np.nansum(poids[:,:thresh,:,:],axis=1)
    poids = (poids1+np.nansum(poids[:,:thresh+1,:,:],axis=1))/2.0

    outmodis = modis*poids
    outmodis[:,latsm>80,:] = np.nan
    
    return outmodis,modmask,latsm,lonsm

aqua,modmaska,latsm,lonsm = open_modis('aqua')
terra,modmask,latsm,lonsm = open_modis('terra')

mods = ['BCC-ESM','EC-Earth','IPSL-CM6','MPI-ESM','MRI-ESM','NorESM']
mmrsf = ['od550ss_AERmon_BCC-ESM1_historical_r1i1p1f1_gn_185001-201412.nc','od550ss_AERmon_EC-Earth3-AerChem_historical_r1i1p1f1_gn_195001-201412.nc','od550ss_AERmon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gr_185001-201412.nc','od550ss_AERmon_MPI-ESM-1-2-HAM_historical_r3i1p1f1_gn_195001-201412.nc','od550ss_AERmon_MRI-ESM2-0_historical_r1i1p1f1_gn_185001-201412.nc','od550ss_AERmon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc']

yy = 10

i=0
for md in mods:
    if i==0:
        init = netCDF4.Dataset('../../../CMIP_DATA/aod/regridded_hr/'+mmrsf[i],'r')
        lons = init.variables['lon'][:]
        lats = init.variables['lat'][:]
        var = init.variables['od550ss'][-yy*12:,:,:]
        var[var>1000] = np.nan
        vv = np.zeros((len(mods),np.shape(var)[0],np.shape(var)[1],np.shape(var)[2]))
        vv[i,:,:,:] = var
    else:
        init = netCDF4.Dataset('../../../CMIP_DATA/aod/regridded_hr/'+mmrsf[i],'r')
        lons = init.variables['lon'][:]
        lats = init.variables['lat'][:]
        vv[i,:,:,:]=init.variables['od550ss'][-yy*12:,:,:]
    i=i+1


vv_t = vv*0
vv_a = vv*0
for hh in range(np.shape(vv_t)[0]):
    vv_t[hh,:] = vv[hh,:]*modmask
    vv_a[hh,:] = vv[hh,:]*modmaska

cmip_t = np.array([np.nanmean(vv_t[:,h::12,:,:],axis=1) for h in range(12)])
arc_t = np.nanmean(cmip_t[:,:,lats>60,:],axis=(2,3))
ant_t =	np.nanmean(cmip_t[:,:,lats<-60,:],axis=(2,3))

cmip_a = np.array([np.nanmean(vv_a[:,h::12,:,:],axis=1) for h in range(12)])
arc_a =	np.nanmean(cmip_a[:,:,lats>60,:],axis=(2,3))
ant_a = np.nanmean(cmip_a[:,:,lats<-60,:],axis=(2,3))

mod_t = np.array([np.nanmean(terra[h::12,:,:],axis=0) for h in range(12)])
mod_arct = np.nanmean(mod_t[:,latsm>60,:],axis=(1,2))
mod_antt = np.nanmean(mod_t[:,latsm<-60,:],axis=(1,2))

mod_a = np.array([np.nanmean(aqua[h::12,:,:],axis=0) for h in range(12)])
mod_arca = np.nanmean(mod_a[:,latsm>60,:],axis=(1,2))
mod_anta = np.nanmean(mod_a[:,latsm<-60,:],axis=(1,2))

arcsdt = np.array([np.nanstd(terra[h::12,latsm>60,:]) for h in range(12)])
arcbct = np.array([np.count_nonzero(~np.isnan(terra[h::12,latsm>60,:])) for h in range(12)])
arcsda = np.array([np.nanstd(aqua[h::12,latsm>60,:]) for h in range(12)])
arcbca = np.array([np.count_nonzero(~np.isnan(aqua[h::12,latsm>60,:])) for h in range(12)])

antsdt = np.array([np.nanstd(terra[h::12,latsm<-60,:]) for h in range(12)])
antbct = np.array([np.count_nonzero(~np.isnan(terra[h::12,latsm<-60,:])) for h in range(12)])
antsda = np.array([np.nanstd(aqua[h::12,latsm<-60,:]) for h in range(12)])
antbca = np.array([np.count_nonzero(~np.isnan(aqua[h::12,latsm<-60,:])) for h in range(12)])

colormap = plt.get_cmap('oranges')
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=[colormap(k) for k in np.linspace(0, 1, 10)])
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
colors=colors[2:-1]

f, axs = plot.subplots(ncols=4,nrows=2,proj={1:'npaeqd',2:'npaeqd',3:'npaeqd',4:'cartesian',5:'spaeqd',6:'spaeqd',7:'spaeqd',8:'cartesian'},sharey=False,sharex=False)

axs[3].bar(range(12),arcbct/np.nanmax(arcbct)*0.15,color='lightgray',alpha=0.5)
axs[7].bar(range(12),antbct/np.nanmax(antbct)*0.15,color='lightgray',alpha=0.5)

for l in range(len(mods)):
    if l==len(mods)-1:
        axs[3].plot(range(12),arc_t[:,l],label='CMIP6',color=colors[l],zorder=9)
        axs[7].plot(range(12),ant_t[:,l],color=colors[l],zorder=9)
    else:
        axs[3].plot(range(12),arc_t[:,l],label='_no_legend_',color=colors[l],zorder=9)
        axs[7].plot(range(12),ant_t[:,l],color=colors[l],zorder=9)
      
axs[3].errorbar(np.arange(-0.1,11.9,1),mod_arct,arcsdt/2,label='Terra',color='black',zorder=10,linewidth=2)
p,q,r = axs[3].errorbar(np.arange(0.1,12.1,1),mod_arca,arcsda/2,label='Aqua',color='grey',zorder=10,linewidth=2)
axs[3].set_xticks(range(12))
axs[3].set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
axs[3].set_ylim(ymin=0.0)
axs[3].set_xlim(-0.2,11.2)
axs[3].legend(ncol=1)
axs[3].set_title('>60$^{\circ}$N')

axs[7].errorbar(np.arange(-0.1,11.9,1),mod_antt,antsdt/2,label='Terra',color='black',zorder=10,linewidth=2)
p,q,r = axs[7].errorbar(np.arange(0.1,12.1,1),mod_anta,antsda/2,label='Aqua',color='grey',zorder=10,linewidth=2)
axs[7].set_xticks(range(12))
axs[7].set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
axs[7].set_ylim(ymin=0.0)
axs[7].set_xlim(-0.2,11.2)
axs[7].set_title('<60$^{\circ}$S')

cm = 'terrain'

vmi,vma = 0,0.16
ttickes = np.linspace(vmi,vma,9)

axs[2].format(coast=True, latlines=10,coastcolor='black',land=True,landcolor='white',title='CMIP6 ensemble',boundinglat=50)
axs[0].format(coast=True, latlines=10,coastcolor='black',land=True,landcolor='white',title='MOD08_M3 (Terra)',boundinglat=50)
axs[1].format(coast=True, latlines=10,coastcolor='black',land=True,landcolor='white',title='MYD08_M3 (Aqua)',boundinglat=50)
axs[4].format(coast=True, latlines=10,coastcolor='black',land=True,landcolor='white',title='',boundinglat=-40)
axs[5].format(coast=True, latlines=10,coastcolor='black',land=True,landcolor='white',title='',boundinglat=-40)
axs[6].format(coast=True, latlines=10,coastcolor='black',land=True,landcolor='white',title='',boundinglat=-40)

m=axs[0].pcolormesh(lonsm,latsm,np.nanmean(terra[:,:,:],axis=0),cmap=cm,vmin=vmi,vmax=vma,alpha=0.9,N=50,extend='max')
m1=axs[4].pcolormesh(lonsm,latsm,np.nanmean(terra[:,:,:],axis=0),cmap=cm,vmin=vmi,vmax=vma,alpha=0.9,N=50,extend='max')
m2=axs[1].pcolormesh(lonsm,latsm,np.nanmean(aqua[:,:,:],axis=0),cmap=cm,vmin=vmi,vmax=vma,alpha=0.9,N=50,extend='max')
m3=axs[5].pcolormesh(lonsm,latsm,np.nanmean(aqua[:,:,:],axis=0),cmap=cm,vmin=vmi,vmax=vma,alpha=0.9,N=50,extend='max')
m__=axs[2].pcolormesh(lons,lats,np.nanmedian(np.nanmean(vv[:,:,:,:],axis=1),axis=(0)),cmap=cm,vmin=vmi,vmax=vma,alpha=0.9,N=50,extend='max')
m_=axs[6].pcolormesh(lons,lats,np.nanmedian(np.nanmean(vv[:,:,:,:],axis=1),axis=(0)),cmap=cm,vmin=vmi,vmax=vma,alpha=0.9,N=50,extend='max')

plt.suptitle('Sea-salt aerosol optical depth at 550nm')
plt.subplots_adjust(hspace=0.05,wspace=0.05)
plt.show()
