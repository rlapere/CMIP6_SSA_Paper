########################
# Remy Lapere (03/11/2022)
# Figure 7-8-AX for Lapere et al., 2022 - CMIP6 SSaer
########################

import warnings
from matplotlib import ticker, cm
import numpy as np
import intake
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
import scipy
import matplotlib as mpl
from matplotlib.tri import Triangulation
import matplotlib.lines as lines
plt.rcParams.update({'font.size': 14})
warnings.filterwarnings('ignore')

df_tot = pd.read_csv('../../../CMIP6/SSA_data_seasonal/total.csv',sep=',',index_col=False)

stats = [['Alert',82.4991,-62.34152],['Bredkalen',63.85,15.33],['Irafoss',64.083,-21.016],['Karasjok',69.46,25.21],['Pallas',68,24.237],['Summit',72.58,-38.479],['Barrow',71.3230,-156.61],['Nord',81.6,-16.67],['Zeppelin',79,12],['Concordia',-75.09,123.33],['Ddu',-66.67,140],['Halley',-75.605,-26.21],['Neumayer',-70.666,-8.266],['Palmer',-64.77,-64.05]]

statof = [['Alert',82.4991,-62.34152],[r'Bredk$\"{a}$len',63.85,15.33],['Irafoss',64.083,-21.016],['Karasjok',69.46,25.21],['Pallas',68,24.237],['Summit',72.58,-38.479],['Utqia$\dot{g}$vik',71.3230,-156.61],['Villum',81.6,-16.67],['Zeppelin',79,12],['Concordia',-75.09,123.33],["Dumont d'U.",-66.67,140],['Halley',-75.605,-26.21],['Neumayer',-70.666,-8.266],['Palmer',-64.77,-64.05]]


lgd = ['BCC','CESM','EC','GISS','HADGEM','IPSL','MIROC','MPI','MRI','NORESM','UKESM1']
lgd1 = ['BCC-ESM','CESM','EC-Earth','GISS','HadGEM','IPSL-CM6','MIROC-ES2L','MPI-ESM','MRI-ESM','NorESM','UKESM','ENS','CAMS','MERRA2']
lgd11 = ['BCC-ESM','CESM','EC-Earth','GISS','HadGEM','IPSL-CM6','MIROC-ES2L','MPI-ESM','MRI-ESM','NorESM','UKESM','Ensemble','CAMS','MERRA2']

models = ['mmrss_AERmon_BCC-ESM1_historical_r1i1p1f1_gn_185001-201412.nc','mmrss_AERmon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','mmrss_AERmon_EC-Earth3-AerChem_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_GISS-E2-1-H_historical_r1i1p3f1_gn_195101-201412.nc','mmrss_AERmon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_195001-201412.nc','mmrss_AERmon_IPSL-CM6A-LR-INCA_historical_r1i1p1f1_gr_185001-201412.nc','mmrss_AERmon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','mmrss_AERmon_MPI-ESM-1-2-HAM_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_MRI-ESM2-0_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_NorESM2-LM_historical_r1i1p1f1_gn_195001-201412.nc','mmrss_AERmon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc']

cor = pd.DataFrame(index=[stat[0] for stat in statof],columns=lgd11,dtype=float)
bias = pd.DataFrame(index=[stat[0] for stat in statof],columns=lgd11,dtype=float)

convnacl = 0.3061 # Seinfeld

def find_point(la,lo,lats_,lons_):
    la_ = np.argmin(np.abs(lats_-la))
    lo_ = np.argmin(np.abs(lons_-lo))
    return la_,lo_


cams_ = netCDF4.Dataset('../../../CAMS_SSA.nc','r')
lons = cams_.variables['longitude'][:]
lats = cams_.variables['latitude'][:]
dens = 1.2922
conc =cams_.variables['aermr01'][:]+cams_.variables['aermr02'][:]+cams_.variables['aermr03'][:]
var = conc*dens*convnacl
camsvar = var/4.3 # reduction factor for dry mass only following ref below
#https://confluence.ecmwf.int/display/CKB/CAMS+global+sea+salt+aerosol+mixing+ratios

merraf = netCDF4.Dataset('../../../MERRA2_SSCONC_MONTHLY.nc','r')
latc = merraf.variables['lat'][:]
lonc = merraf.variables['lon'][:]
merr = merraf.variables['SSSMASS'][:]*convnacl

f,axs = pplt.subplots(ncols=5,nrows=3,sharey=0)
f1,axs1 = pplt.subplots(ncols=5,nrows=3,sharey=0)

ind=1

for st in stats:
    obstat = df_tot[st[0]]
    sdstat = df_tot[st[0]+'_sd']
    i=0
    vv = np.zeros((len(models),12))
    vv1 = np.zeros((len(models),12))
    for md in models:
        df_cmip = pd.read_csv('EXTCMIP/'+st[0]+'_'+lgd1[i]+'.csv')
        cmip = df_cmip['ssa'].values*dens*convnacl
        if np.isin(md,['IPSL-CM6']):
            cmip = cmip/4.3
        vv[i,:] = cmip
        biass = np.mean(cmip)/np.mean(obstat)
        vv1[i,:] = cmip/biass
        bias.iloc[ind-1, i] = (np.mean(vv[i,:])-np.mean(obstat))/np.mean(obstat)

        alpha = 0.05
        z_critical = scipy.stats.norm.ppf(1 - alpha)
        r=pearsonr(obstat,vv1[i,:])[0]
        z_prime = 0.5*np.log((1+r)/(1-r))
        n=np.shape(vv1[i,:])[0]
        se = 1./np.sqrt(n-3)
        z_low = z_prime-z_critical*se
        z_up = z_prime+z_critical*se
        CI_1 = (np.exp(2*z_low)-1)/(np.exp(2*z_low)+1)
        CI_0 = (np.exp(2*z_up)-1)/(np.exp(2*z_up)+1)
        cor.iloc[ind-1, i] = r
        i=i+1
        
    ens = np.mean(vv,axis=0)
    med = np.median(vv,axis=0)
    ens1 = np.mean(vv1,axis=0)

    bias.iloc[ind-1, i] = (np.mean(ens)-np.mean(obstat))/np.mean(obstat)

    alpha = 0.05
    z_critical = scipy.stats.norm.ppf(1 - alpha)
    r=pearsonr(obstat,ens1)[0]
    z_prime = 0.5*np.log((1+r)/(1-r))
    n=np.shape(ens1)[0]
    se = 1./np.sqrt(n-3)
    z_low = z_prime-z_critical*se
    z_up = z_prime+z_critical*se
    CI_1 = (np.exp(2*z_low)-1)/(np.exp(2*z_low)+1)
    CI_0 = (np.exp(2*z_up)-1)/(np.exp(2*z_up)+1)
    cor.iloc[ind-1, i] = r
    cortex = str(round(r,2))+' ['+str(round(CI_1,2))+','+str(round(CI_0,2))+']'
    i=i+1
    
    if st[2]<0:
        stloc = st[2]+360
    else:
        stloc = st[2]
    x,y = find_point(st[1],stloc,lats,lons)    
    camsloc = camsvar[:,x,y]
    cams = np.array([np.mean(camsloc[i::12]) for i in range(12)])*10**9    
    camssd = np.array([np.std(camsloc[i::12]) for i in range(12)])*10**9
    cor.iloc[ind-1, i] = pearsonr(obstat,cams)[0]
    bias.iloc[ind-1, i] = (np.mean(cams)-np.mean(obstat))/np.mean(obstat)
    i=i+1

    x,y = find_point(st[1],st[2],latc,lonc)
    merloc = merr[:,x,y]
    merra = np.array([np.mean(merloc[i::12]) for i in range(12)])*10**9
    merrasd = np.array([np.std(merloc[i::12]) for i in range(12)])*10**9
    cor.iloc[ind-1, i] = pearsonr(obstat,merra)[0]
    bias.iloc[ind-1, i] = (np.mean(merra)-np.mean(obstat))/np.mean(obstat)

    if ind>9:
        indloc=ind+1
    else:
        indloc=ind
        
    if ind==1:
        labs=['CMIP6 mean','CMIP6 median','CAMS','MERRA2','OBS']
    else:
        labs=['_no_legend_','_no_legend_','_no_legend_','_no_legend_','_no_legend_']

    colormap = plt.get_cmap('blues6')

    cole = colormap(0.85)
    coli = colormap(0.6)
    cole = colormap(0.7)
    colr = 'chocolate'
    axs[indloc-1].errorbar(range(12),obstat,sdstat/2,zorder=10,color='black',label=labs[4])
    axs[indloc-1].plot(vv.T,color=coli,zorder=0,alpha=0.5,linewidth=0.5)
    axs[indloc-1].plot(ens,color=cole,zorder=3,label=labs[0],linewidth=2.5)
    axs[indloc-1].plot(med,color=cole,zorder=3,label=labs[1],linewidth=2.5,linestyle='--')
    axs[indloc-1].plot(range(12),cams,color=colr,zorder=1,label=labs[2],marker='o',markersize=4,alpha=0.5)
    axs[indloc-1].plot(range(12),merra,color=colr,zorder=1,label=labs[3],marker='^',markersize=4,alpha=0.5)
    axs[indloc-1].set_yscale('log')
    axs[indloc-1].set_ylim(ymin=0)
    axs[indloc-1].set_title(statof[ind-1][0])
    
    if ind<=9:
        axs[indloc-1].set_xticks([])
    else:
        axs[indloc-1].set_xticks(range(12))
        axs[indloc-1].set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])

    if np.mod(indloc-1,5)==0:
        axs[indloc-1].set_ylabel('\u03BC'+'g/m$^3$')

    axs1[indloc-1].errorbar(range(12),obstat,sdstat/2,zorder=10,color='black',label=labs[4])
    axs1[indloc-1].plot(range(12),vv1.T,color=coli,zorder=0,alpha=0.5,linewidth=0.5)
    axs1[indloc-1].plot(range(12),ens/np.mean(ens)*np.mean(obstat),color=cole,zorder=1,label=labs[0],linewidth=2.5)
    axs1[indloc-1].set_title(statof[ind-1][0])
    axs1[indloc-1].text(0.4,0.9,s=cortex,size=9,bbox={'edgecolor':'black','facecolor':'white'},zorder=11,transform=axs1[indloc-1].transAxes)
    axs1[indloc-1].set_ylim(ymin=0)

    if ind<=9:
        axs1[indloc-1].set_xticks([])
    else:
        axs1[indloc-1].set_xticks(range(12))
        axs1[indloc-1].set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])

    if np.mod(indloc-1,5)==0:
        axs1[indloc-1].set_ylabel('\u03BC'+'g/m$^3$')

    ind = ind+1

f.legend(loc='bottom',fontsize=10,ncol=5)
axs.format(suptitle='Modeled and observed mass concentration')

f1.legend(loc='bottom',fontsize=10,ncol=5)
axs1.format(suptitle='Modeled and observed mass concentration')

plt.show()

statof = ['Alert','Bredkalen','Irafoss','Karasjok','Pallas','Summit','Utqiagvik','Villum','Zeppelin','Concordia',"Dumont d'U.",'Halley','Neumayer','Palmer']
gd11 = ['BCC-ESM','CESM','EC-Earth','GISS','HadGEM','IPSL-CM6','MIROC-ES2L','MPI-ESM','MRI-ESM','NorESM','UKESM','Ensemble','CAMS','MERRA2']
cor['ensmed'] = np.zeros(len(cor.index))

for kk in range(len(cor.index)):
    cor['ensmed'][kk] = len(np.where(cor.iloc[kk]>0)[0])

bias['ensco'] = cor.ensmed.values
bias['ens'] = cor.Ensemble.values
cor = cor.sort_values(['ensmed'])
bias = bias.sort_values(['ensco'])
bias = bias.drop(labels=["ensco",'ens'],axis=1)
cor = cor.drop(labels='ensmed',axis=1)

ax=sns.heatmap(cor.to_numpy().astype(float),annot=bias.to_numpy().astype(float)*100,cmap='ColdHot',vmin=-1,vmax=1,label='Correlation',fmt='.0f',cbar_kws={'label':'Correlation'},linewidths=0.5,linecolor='k')

ax.set_xticks(np.arange(0.5,len(gd11)+0.5,1))
ax.set_xticklabels(gd11,rotation=90)
ax.set_yticks(np.arange(0.5,len(cor.Ensemble.values)+0.5,1))
ax.set_yticklabels(cor.index.values,rotation=0,color='k')

col=np.array(['k','k','k','k','k','k','teal','teal','k','k','teal','teal','teal','k'])
col=np.flip(col)

i=0
for t in ax.yaxis.get_ticklabels():
    t.set_color(col[i])
    i=i+1

plt.show()
