########################
# Remy Lapere (03/11/2022)
# Figure 10 for Lapere et al., 2022 - CMIP6 SSaer
########################

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,BoundaryNorm,LinearSegmentedColormap
import netCDF4
import numpy as np
import proplot as plot
import pandas as pd
import cartopy.crs as ccrs
import cartopy
import matplotlib as mpl
import scipy
from scipy.stats import wilcoxon,norm
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


fld = '/Users/rlapere/Desktop/CMIP_SSA/CMIP_DATA/aod/'

mods = ['a) IPSL-CM6','b) NorESM','c) UKESM']

ctrl = ['od550aer_AERmon_IPSL-CM6A-LR-INCA_piClim-control_r1i1p1f1_gr_185001-187912.nc','od550aer_AERmon_NorESM2-LM_piClim-control_r1i1p1f1_gn_000101-003012.nc','od550aer_AERmon_UKESM1-0-LL_piClim-control_r1i1p1f4_gn_185001-189412.nc']
xss = ['od550aer_AERmon_IPSL-CM6A-LR-INCA_piClim-2xss_r1i1p1f1_gr_185001-187912.nc','od550aer_AERmon_NorESM2-LM_piClim-2xss_r1i1p1f1_gn_000101-003012.nc','od550aer_AERmon_UKESM1-0-LL_piClim-2xss_r1i1p1f4_gn_185001-189412.nc']

fig, axs = plot.subplots(ncols=3,nrows=4,proj={1:'npaeqd',2:'npaeqd',3:'npaeqd',4:'spaeqd',5:'spaeqd',6:'spaeqd',7:'npaeqd',8:'npaeqd',9:'npaeqd',10:'spaeqd',11:'spaeqd',12:'spaeqd'},sharey=False,sharex=False)
axs.format(coast=True, latlines=10, coastcolor='black', suptitle='Change in radiation with doubled sea-salt aerosol emissions')

k=0
for md in mods:
    print(md)
    rtmt = netCDF4.Dataset(fld+ctrl[k],'r')
    rtmt2 = netCDF4.Dataset(fld+xss[k],'r')
    rt2 = rtmt2.variables['od550aer'][-30*12:,:,:]
    rt = rtmt.variables['od550aer'][-30*12:,:,:]

    seas = range(len(rtmt2.variables['od550aer'][:,0,0]))
    hiv = np.where((np.mod(seas,12)==0) | (np.mod(seas,12)==1) | (np.mod(seas,12)==11))[0]
    ete = np.where((np.mod(seas,12)==5) | (np.mod(seas,12)==6) | (np.mod(seas,12)==7))[0]
    spr = np.where((np.mod(seas,12)==2) | (np.mod(seas,12)==3) | (np.mod(seas,12)==4))[0]
    fal = np.where((np.mod(seas,12)==8) | (np.mod(seas,12)==9) | (np.mod(seas,12)==10))[0]
    
    lats = rtmt.variables['lat'][:]
    lons = rtmt.variables['lon'][:]

    rt_djf = rt[hiv,:,:]
    rt2_djf = rt2[hiv,:,:]
    rt_jja = rt[ete,:,:]
    rt2_jja = rt2[ete,:,:]
    
    indn = np.where(lats>=60)[0]
    inds = np.where(lats<=-60)[0]

    gwn = gridw(lats[indn],lons)
    gws = gridw(lats[inds],lons)
    
    difarc_djf = np.mean(rt2_djf[:,indn,:]-rt_djf[:,indn,:],axis=0).flatten()*gwn.flatten()
    difant_djf = np.mean(rt2_djf[:,inds,:]-rt_djf[:,inds,:],axis=0).flatten()*gws.flatten()
    difarc_jja = np.mean(rt2_jja[:,indn,:]-rt_jja[:,indn,:],axis=0).flatten()*gwn.flatten()
    difant_jja = np.mean(rt2_jja[:,inds,:]-rt_jja[:,inds,:],axis=0).flatten()*gws.flatten()

    alpha = 0.05
    z = norm.ppf(1-alpha)
    ciarc_djf = z*np.std(difarc_djf)/np.sqrt(len(difarc_djf))
    ciant_djf = z*np.std(difant_djf)/np.sqrt(len(difant_djf))
    ciarc_jja = z*np.std(difarc_jja)/np.sqrt(len(difarc_jja))
    ciant_jja = z*np.std(difant_jja)/np.sqrt(len(difant_jja))
    
    wilcodjf = np.zeros((len(lats),len(lons)))*np.nan
    wilcojja = np.zeros((len(lats),len(lons)))*np.nan
    i=0
    for la in lats:
        if (la>45) or (la<-50):
            j=0
            for lo in lons:
                w,wilcodjf[i,j] = wilcoxon(rt2_djf[:,i,j]-rt_djf[:,i,j])
                w,wilcojja[i,j] = wilcoxon(rt2_jja[:,i,j]-rt_jja[:,i,j])
                j=j+1
        i=i+1

    wilcodjf[wilcodjf>0.1]=np.nan
    wilcodjf[~np.isnan(wilcodjf)] = 1

    wilcojja[wilcojja>0.1]=np.nan
    wilcojja[~np.isnan(wilcojja)] = 1

    diff_djf = np.mean(rt2_djf-rt_djf,axis=0)
    diff_jja = np.mean(rt2_jja-rt_jja,axis=0)

    x,y = np.meshgrid(lons,lats)
    x[wilcodjf==1] = np.nan
    x1,y1 = np.meshgrid(lons,lats)
    x1[wilcojja==1] = np.nan

    m3=axs[k].pcolormesh(lons,lats,diff_djf,cmap='lajolla',vmin=0,vmax=0.25,N=10)
    axs[k].plot(x,y,linewidth=0,marker='x',markersize=1,markerfacecolor=(0.25,0.25,0.25,0.75),markeredgecolor=(0.25,0.25,0.25,0.75))

    m=axs[k+3].pcolormesh(lons,lats,diff_jja,cmap='lajolla',vmin=0,vmax=0.25,N=10)
    axs[k+3].plot(x1,y1,linewidth=0,marker='x',markersize=1,markerfacecolor=(0.25,0.25,0.25,0.75),markeredgecolor=(0.25,0.25,0.25,0.75))

    m1=axs[k+6].pcolormesh(lons,lats,diff_jja,cmap='lajolla',vmin=0,vmax=0.25,N=10)
    axs[k+6].plot(x1,y1,linewidth=0,marker='x',markersize=1,markerfacecolor=(0.25,0.25,0.25,0.75),markeredgecolor=(0.25,0.25,0.25,0.75))

    m2=axs[k+9].pcolormesh(lons,lats,diff_djf,cmap='lajolla',vmin=0,vmax=0.25,N=10)
    axs[k+9].plot(x,y,linewidth=0,marker='x',markersize=1,markerfacecolor=(0.25,0.25,0.25,0.75),markeredgecolor=(0.25,0.25,0.25,0.75))

    axs[k].set_title(md)

    k=k+1

axs[0].format(boundinglat=50)
axs[1].format(boundinglat=50)
axs[2].format(boundinglat=50)
axs[3].format(boundinglat=-50)
axs[4].format(boundinglat=-50)
axs[5].format(boundinglat=-50)
axs[6].format(boundinglat=50)
axs[7].format(boundinglat=50)
axs[8].format(boundinglat=50)
axs[9].format(boundinglat=-50)
axs[10].format(boundinglat=-50)
axs[11].format(boundinglat=-50)
plt.show()
