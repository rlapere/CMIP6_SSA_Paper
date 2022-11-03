########################
# Remy Lapere (03/11/2022)
# Figure 1 for Lapere et al., 2022 - CMIP6 SSaer
########################

import matplotlib.pyplot as plt
import proplot as pplt
import netCDF4
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy
from matplotlib.colors import LinearSegmentedColormap

colsta = np.repeat('white',9)
colsta1=np.repeat('white',5)

stats1 = [['Neumayer',-70.66,-8.26],['Halley',-75.60,-26.21],['Concordia',-75.09,123.33],["Dumont d'U.",-66.67,140],['Palmer',-64.77,-64.05]]
stats = [['Irafoss',64.08,-21.01],['Bredkalen',63.85,15.33],['Karasjok',69.46,25.21],['Pallas',68.00,24.23],['Zeppelin',79,12],['Villum',81.60,-16.67],['Alert',82.49,-62.34],['Utqia$\dot{g}$vik',71.32,-156.61],['Summit',72.58,-38.48]]

stations = pd.DataFrame({'Station':np.append([st[0] for st in stats],[st[0] for st in stats1]),
                         'Lat':np.append([st[1] for st in stats],[st[1] for st in stats1]),
                         'Lon':np.append([st[2] for st in stats],[st[2] for st in stats1]),
                         'Period':['2006-2020','2009-2020','2001-2010','2003-2019','1993-2020','2000-2021','1980-1990','1984-2001','2003-2006','1983-2007','1991-1993','1990-1998','1991-2013','1990-1996'],
                         'Source':['EBAS','EBAS','EBAS','EBAS','EBAS','EBAS','Yang et al., 2019','Yang et al., 2019','Yang et al., 2019','Yang et al., 2019','Yang et al., 2019','Yang et al., 2019','Legrand et al., 2016','Yang et al., 2019']})

seaice = netCDF4.Dataset('/Users/rlapere/Downloads/ERA5_siconc_hist.nc')
si = seaice.variables['siconc'][-20:,0,:,:]
lats = seaice.variables['latitude'][:]
lons = seaice.variables['longitude'][:]

si[si<0] = 1
simin = np.nanmean(si[1::2,:,:],axis=0)
simax = np.nanmean(si[0::2,:,:],axis=0)
simin[simin>0.5] = 100
simin[simin<=0.5] = np.nan
simax[simax>0.5] = 100
simax[simax<=0.5] = np.nan

cmap = LinearSegmentedColormap.from_list('testCmap', [(1,1,1,1),(181/255.,211/255.,231/255.,1)], N=2)
cmap1 = LinearSegmentedColormap.from_list('testCmap', [(1,1,1,1),(134/255.,161/255.,183/255.,1)], N=2)

ax1 = plt.subplot(1,2,1,projection=ccrs.NorthPolarStereo())                                                                                                                 
ax1.set_extent([-180,180,58,90], ccrs.PlateCarree())                                                                                                                                                       
ax1.add_feature(cartopy.feature.LAND,edgecolor='k',zorder=4,facecolor=(0.5,0.5,0.5,1))

c=ax1.pcolormesh(lons,lats,simin,zorder=3,transform=ccrs.PlateCarree(),vmin=10,vmax=100,alpha=0.5,cmap=cmap1,shading='gouraud')
c.cmap.set_under((1,1,1,0))
c0=ax1.pcolormesh(lons,lats,simax,zorder=2,transform=ccrs.PlateCarree(),vmin=10,vmax=100,alpha=0.7,cmap=cmap,shading='gouraud')
c0.cmap.set_under((1,1,1,0))

ax1.plot(np.arange(-180,181,1),np.repeat(60,len(np.arange(-180,181,1))),linestyle='--',color='k',transform=ccrs.PlateCarree(),zorder=5,linewidth=0.5)

k=0                                                                                                                                                                                                       
for stat in stats:                                                                                                                                                                                
    ax1.text(stat[2],stat[1],stat[0][:2],transform=ccrs.PlateCarree(),bbox=dict(facecolor=colsta[k]),zorder=6)
    k=k+1

statar = stations[stations.Lat>0]
table = ax1.table(cellText=statar.values, colLabels=statar.columns,colWidths=[0.2,0.16,0.16,0.2,0.28],cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(10)

ax = plt.subplot(1,2,2,projection=ccrs.SouthPolarStereo())                                                                                                                                                 
ax.set_extent([-180,180,-58,-90], ccrs.PlateCarree())                                                                                                                                                      
ax.add_feature(cartopy.feature.LAND,edgecolor='k',zorder=4,facecolor=(0.5,0.5,0.5,1))                                                                                                                      
c1=ax.pcolormesh(lons,lats,simin,zorder=2,transform=ccrs.PlateCarree(),vmin=10,vmax=100,alpha=0.7,cmap=cmap,shading='gouraud')
c1.cmap.set_under((1,1,1,0))
c10=ax.pcolormesh(lons,lats,simax,zorder=3,transform=ccrs.PlateCarree(),vmin=10,vmax=100,alpha=0.5,cmap=cmap1,shading='gouraud')
c10.cmap.set_under((1,1,1,0))

ax.plot(np.arange(-180,181,1),np.repeat(-60,len(np.arange(-180,181,1))),linestyle='--',color='k',transform=ccrs.PlateCarree(),zorder=5,linewidth=0.5)


k=0
for stat in stats1:                                                                                                                                                                                        
    ax.text(stat[2],stat[1],stat[0][:2],transform=ccrs.PlateCarree(),bbox=dict(facecolor=colsta1[k]),zorder=6)                                                                                                      
    k=k+1                                                                                                                                                                                                  

statan = stations[stations.Lat<0]
table = ax.table(cellText=statan.values, colLabels=statan.columns,colWidths=[0.2,0.14,0.14,0.2,0.32],cellLoc='center')

table.auto_set_font_size(False)
table.set_fontsize(10)

plt.subplots_adjust(left=0.02,
                    bottom=0.35,
                    right=0.99,
                    top=0.98,
                    wspace=0.01,
                    hspace=0.01)

plt.show()
