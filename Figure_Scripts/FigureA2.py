########################
# Remy Lapere (03/11/2022)
# Figure A2 for Lapere et al., 2022 - CMIP6 SSaer
########################

import warnings
from matplotlib import ticker, cm
import numpy as np
import intake
import xarray as xr 
import proplot as pplt 
import matplotlib.pyplot as plt
from intake import open_catalog
import cartopy.crs as ccrs
import cartopy
import pandas as pd
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy.stats import pearsonr
from matplotlib.colors import LogNorm
import netCDF4
plt.rcParams.update({'font.size': 12})
warnings.filterwarnings('ignore')

df_tot = pd.read_csv('../../../CMIP6/SSA_data_seasonal/total.csv',sep=',',index_col=False)
obstat = df_tot['Alert']
sdstat = df_tot['Alert_sd']

mods = ['CESM','EC-Earth','GISS','HadGEM','IPSL-CM6','MIROC-ES2L','MPI-ESM','NorESM','UKESM']

col=['red','green','blue','yellow','black','gray','orange','pink','yellow']
col=['cadetblue','cadetblue','cadetblue','pink','cadetblue','cadetblue','cadetblue','cadetblue','crimson']
col=['gray','gray','gray','pink','gray','gray','gray','gray','crimson']

f,axs = pplt.subplots(ncols=1,nrows=3,sharey=0)
ax0 = axs[0]
ax = axs[1]
ax1 = axs[2]

ax0.set_title('Normalized SSaer surface mass concentration')
ax.set_title('Normalized SSaer mass emission rate')
ax1.set_title('Planetary boundary layer height')

mini = np.zeros((3,12))+1000
maxi = np.zeros((3,12))

i=0
vv = np.zeros((len(mods),12))
for md in mods:
    df_emi = pd.read_csv('EXTCMIP/Alert_'+md+'_emis.csv')
    df_pbl = pd.read_csv('EXTCMIP/Alert_'+md+'_bldep.csv')
    df_ = pd.read_csv('EXTCMIP/Alert_'+md+'.csv')
    if i>0:
        lab='_no_legend_'
        lw=1
        ap=0.6
    if i==0:
        lab='Other models'
        lw=1
        ap=0.6
    if md=='HadGEM':
        lab='HadGEM'
        lw=1.75
        ap=1
    if md=='UKESM':
        lab='UKESM'
        lw=1.75
        ap=1

    ax.plot(range(12),df_emi.ssa.values/np.max(df_emi.ssa.values),color=col[i],linewidth=lw,alpha=ap)
    ax1.plot(range(12),df_pbl.ssa.values,color=col[i],linewidth=lw,alpha=ap)
    ax0.plot(range(12),df_.ssa.values/np.max(df_.ssa.values),color=col[i],label=lab,linewidth=lw,alpha=ap)
    
    if(~np.isin(md,['HadGEM','UKESM'])):
        mini[0,:] = np.fmin(mini[0,:],df_.ssa.values/np.max(df_.ssa.values))
        maxi[0,:] =	np.fmax(maxi[0,:],df_.ssa.values/np.max(df_.ssa.values))
        mini[1,:] = np.fmin(mini[1,:],df_emi.ssa.values/np.max(df_emi.ssa.values))
        maxi[1,:] =	np.fmax(maxi[1,:],df_emi.ssa.values/np.max(df_emi.ssa.values))
        mini[2,:] =	np.fmin(mini[2,:],df_pbl.ssa.values)
        maxi[2,:] = np.fmax(maxi[2,:],df_pbl.ssa.values)

    i=i+1
        
ax.fill_between(range(12),mini[1,:],maxi[1,:],zorder=0,color=(0.5,0.5,0.5,0.25))
ax0.fill_between(range(12),mini[0,:],maxi[0,:],zorder=0,color=(0.5,0.5,0.5,0.25))
ax1.fill_between(range(12),mini[2,:],maxi[2,:],zorder=0,color=(0.5,0.5,0.5,0.25))

f.legend(loc='bottom',ncol=3)

ax.set_xticks([])
ax0.set_xticks([])
ax1.set_xticks(range(12))
ax1.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.set_xlim(0,11)
ax0.set_xlim(0,11)
ax1.set_xlim(0,11)
ax.set_ylim(ymin=0)
ax0.set_ylim(ymin=0)
ax0.set_ylabel('Normalized mass concentration')
ax.set_ylabel('Normalized mass emission')
ax1.set_ylabel('PBLH (m)')
plt.show()
