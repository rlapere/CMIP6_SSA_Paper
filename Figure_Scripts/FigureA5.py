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
import seaborn as sns
from haversine import haversine, Unit
mpl.rcParams.update({'font.size': 12})

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

fld = '/Users/rlapere/Desktop/CMIP_SSA/CMIP_DATA/tas/'

mods = ['NorESM','IPSL-CM6','UKESM']
ctrl = ['tas_Amon_NorESM2-LM_piClim-control_r1i1p1f1_gn_000101-003012.nc','tas_Amon_IPSL-CM6A-LR-INCA_piClim-control_r1i1p1f1_gr_185001-187912.nc','tas_Amon_UKESM1-0-LL_piClim-control_r1i1p1f4_gn_185001-189412.nc']
xss = ['tas_Amon_NorESM2-LM_piClim-2xss_r1i1p1f1_gn_000101-003012.nc','tas_Amon_IPSL-CM6A-LR-INCA_piClim-2xss_r1i1p1f1_gr_185001-187912.nc','tas_Amon_UKESM1-0-LL_piClim-2xss_r1i1p1f4_gn_185001-189412.nc']
xdust = ['tas_Amon_NorESM2-LM_piClim-2xdust_r1i1p1f1_gn_000101-003012.nc','tas_Amon_IPSL-CM6A-LR-INCA_piClim-2xdust_r1i1p1f1_gr_185001-187912.nc','tas_Amon_UKESM1-0-LL_piClim-2xdust_r1i1p1f4_gn_185001-189412.nc']

modo=[]
seas1=[]
toto=[]

modod=[]
seas1d=[]
totod=[]

k=0
for md in mods:
    print(md)
    rtmt = netCDF4.Dataset(fld+ctrl[k],'r')
    rtmt2 = netCDF4.Dataset(fld+xss[k],'r')
    rtmdu = netCDF4.Dataset(fld+xdust[k],'r')
    rt2 = rtmt2.variables['tas'][:,:,:]
    rt = rtmt.variables['tas'][:,:,:]
    rtd = rtmdu.variables['tas'][:,:,:]
    
    seas = range(len(rt[:,0,0]))
    hiv = np.where((np.mod(seas,12)==0) | (np.mod(seas,12)==1) | (np.mod(seas,12)==11))[0]
    ete = np.where((np.mod(seas,12)==5) | (np.mod(seas,12)==6) | (np.mod(seas,12)==7))[0]
    spr = np.where((np.mod(seas,12)==2) | (np.mod(seas,12)==3) | (np.mod(seas,12)==4))[0]
    fal = np.where((np.mod(seas,12)==8) | (np.mod(seas,12)==9) | (np.mod(seas,12)==10))[0]

    lats = rtmt.variables['lat'][:]
    lons = rtmt.variables['lon'][:]

    indn = np.where(lats>=60)[0]
    inds = np.where(lats<=-60)[0]
    gwn = gridw(lats[indn],lons)
    gws = gridw(lats[inds],lons)
    
    rt2n = rt2[:,indn,:]
    rt2s = rt2[:,inds,:]
    rtn = rt[:,indn,:]
    rts = rt[:,inds,:]
    rt2dn = rtd[:,indn,:]
    rt2ds = rtd[:,inds,:]
    
    difarete = np.mean((rt2n[ete,:,:]-rtn[ete,:,:]),axis=0).flatten()*gwn.flatten()
    difanete = np.mean((rt2s[hiv,:,:]-rts[hiv,:,:]),axis=0).flatten()*gws.flatten()
    difarhiv = np.mean((rt2n[hiv,:,:]-rtn[hiv,:,:]),axis=0).flatten()*gwn.flatten()
    difanhiv = np.mean((rt2s[ete,:,:]-rts[ete,:,:]),axis=0).flatten()*gws.flatten()

    difareted = np.mean((rt2dn[ete,:,:]-rtn[ete,:,:]),axis=0).flatten()*gwn.flatten()
    difaneted = np.mean((rt2ds[hiv,:,:]-rts[hiv,:,:]),axis=0).flatten()*gws.flatten()
    difarhivd = np.mean((rt2dn[hiv,:,:]-rtn[hiv,:,:]),axis=0).flatten()*gwn.flatten()
    difanhivd = np.mean((rt2ds[ete,:,:]-rts[ete,:,:]),axis=0).flatten()*gws.flatten()
    
    difarete = difarete[~np.isnan(difarete)]
    difanete = difanete[~np.isnan(difanete)]
    difarhiv = difarhiv[~np.isnan(difarhiv)]
    difanhiv = difanhiv[~np.isnan(difanhiv)]

    difareted = difareted[~np.isnan(difareted)]
    difaneted = difaneted[~np.isnan(difaneted)]
    difarhivd = difarhivd[~np.isnan(difarhivd)]
    difanhivd = difanhivd[~np.isnan(difanhivd)]
    
    difete = np.append(difarete,difanete)
    difhiv = np.append(difarhiv,difanhiv)
    toloc = np.append(difete,difhiv)
    toto = np.append(toto,toloc)

    difeted = np.append(difareted,difaneted)
    difhivd = np.append(difarhivd,difanhivd)
    tolocd = np.append(difeted,difhivd)
    totod = np.append(totod,tolocd)
    
    arce = np.repeat(md+'arc',len(difarete))
    ante = np.repeat(md+'ant',len(difanete))
    arch = np.repeat(md+'arc',len(difarhiv))
    anth = np.repeat(md+'ant',len(difanhiv))
    un = np.append(arce,ante)
    de = np.append(arch,anth)
    aglo = np.append(un,de)
    modo = np.append(modo,aglo)
    seaete = np.repeat('Summer',len(difarete))
    seahiv = np.repeat('Winter',len(difarhiv))
    seaeteq = np.repeat('Summer',len(difanete))
    seahivq = np.repeat('Winter',len(difanhiv))
    seaagg = np.append(seahiv,seahivq)
    seaagg1 = np.append(seaete,seaeteq)
    ss = np.append(seaagg1,seaagg)
    seas1 = np.append(seas1,ss)

    arced = np.repeat(md+'arc',len(difareted))
    anted = np.repeat(md+'ant',len(difaneted))
    archd = np.repeat(md+'arc',len(difarhivd))
    anthd = np.repeat(md+'ant',len(difanhivd))
    und = np.append(arced,anted)
    ded = np.append(archd,anthd)
    aglod = np.append(und,ded)
    modod = np.append(modod,aglod)
    seaeted = np.repeat('Summer',len(difareted))
    seahivd = np.repeat('Winter',len(difarhivd))
    seaeteqd = np.repeat('Summer',len(difaneted))
    seahivqd = np.repeat('Winter',len(difanhivd))
    seaaggd = np.append(seahivd,seahivqd)
    seaagg1d = np.append(seaeted,seaeteqd)
    ssd = np.append(seaagg1d,seaaggd)
    seas1d = np.append(seas1d,ssd)
    
    k=k+1

agreg = pd.DataFrame({'Change in air surface temperature ('+u'\N{DEGREE SIGN}'+'C)':toto,'Model':modo,'Season':seas1})
agreg['pole'] = np.array([str(agreg['Model'].values[i])[-3:] for i in range(len(seas1))])
agreg.pole[agreg.pole=='arc']='Arctic (>60N)'
agreg.pole[agreg.pole=='ant']='Antarctic (<60S)'
agreg['modmod'] = np.array([str(agreg['Model'].values[i])[:-3] for i in range(len(seas1))])
agreg = agreg.sort_values(by=['Season','modmod'])


agregd = pd.DataFrame({'Change in air surface temperature ('+u'\N{DEGREE SIGN}'+'C)':totod,'Model':modod,'Season':seas1d})
agregd['pole'] = np.array([str(agregd['Model'].values[i])[-3:] for i in range(len(seas1d))])
agregd.pole[agregd.pole=='arc']='Arctic (>60N)'
agregd.pole[agregd.pole=='ant']='Antarctic (<60S)'
agregd['modmod'] = np.array([str(agregd['Model'].values[i])[:-3] for i in range(len(seas1d))])
agregd = agregd.sort_values(by=['Season','modmod'])

c = (181/255.,211/255.,231/255.,1)
c1 =(134/255.,161/255.,183/255.,1)


def plt_violin(df_in,ax):
    i=0
    for mm in mods:
        dfloc = df_in[df_in.modmod==mm]
        sns.violinplot(data=dfloc,x="pole",y='Change in air surface temperature ('+u'\N{DEGREE SIGN}'+'C)',hue="Season",cut=0,split=True,scale='area',palette=[c,c1],inner=None,width=0.9,linewidth=0.25)

        difetel = dfloc[dfloc.Season=='Summer']
        difarete1 = difetel[difetel.pole=='Arctic (>60N)']
        difarete1 = difarete1['Change in air surface temperature ('+u'\N{DEGREE SIGN}'+'C)'].values
        difanete1 = difetel[difetel.pole=='Antarctic (<60S)']
        difanete1 = difanete1['Change in air surface temperature ('+u'\N{DEGREE SIGN}'+'C)'].values
        
        difhivl = dfloc[dfloc.Season=='Winter']
        difarhiv1 = difhivl[difhivl.pole=='Arctic (>60N)']
        difarhiv1 = difarhiv1['Change in air surface temperature ('+u'\N{DEGREE SIGN}'+'C)'].values
        difanhiv1 = difhivl[difhivl.pole=='Antarctic (<60S)']
        difanhiv1 = difanhiv1['Change in air surface temperature ('+u'\N{DEGREE SIGN}'+'C)'].values
        
        alpha = 0.05
        z = norm.ppf(1-alpha)
        ciarc1 = z*np.std(difarete1)/np.sqrt(len(difarete1))
        ciant1 = z*np.std(difanete1)/np.sqrt(len(difanete1))
        ciarch1 = z*np.std(difarhiv1)/np.sqrt(len(difarhiv1))
        cianth1 = z*np.std(difanhiv1)/np.sqrt(len(difanhiv1))

        plt.table([mm,str(round(np.mean(difarete1),3))+u'\N{DEGREE SIGN}'+'C +/- '+str(round(ciarc1,3)),str(round(np.mean(difarhiv1),3))+u'\N{DEGREE SIGN}'+'C +/- '+str(round(ciarch1,3)),mm,str(round(np.mean(difanete1),3))+u'\N{DEGREE SIGN}'+'C +/- '+str(round(ciant1,3)),str(round(np.mean(difanhiv1),3))+u'\N{DEGREE SIGN}'+'C +/- '+str(round(cianth1,3))])
        
        plt.text(-0.48,0.8-0.1*i,mm,size=10)
        plt.text(-0.3,0.8-0.1*i,str(round(np.mean(difarete1),3))+u'\N{DEGREE SIGN}'+'C +/- '+str(round(ciarc1,3)),size=10)
        plt.text(0.1,0.8-0.1*i,str(round(np.mean(difarhiv1),3))+u'\N{DEGREE SIGN}'+'C +/- '+str(round(ciarch1,3)),size=10)
        plt.text(1-0.3,0.8-0.1*i,str(round(np.mean(difanete1),3))+u'\N{DEGREE SIGN}'+'C +/- '+str(round(ciant1,3)),size=10)
        plt.text(1+0.1,0.8-0.1*i,str(round(np.mean(difanhiv1),3))+u'\N{DEGREE SIGN}'+'C +/- '+str(round(cianth1,3)),size=10)
        
        i=i+1

    for violin in ax.collections:
        violin.set_alpha(0.5)
    plt.ylim(-1,1)
    plt.grid(axis='y',which='both',color='gray')
    plt.xlabel('')
    plt.hlines(0,-1,2,zorder=10,linewidths=1,colors='k')
    plt.legend([],[], frameon=False)
    plt.xlim(-0.5,1.5)


ax=plt.subplot(121)
plt_violin(agreg,ax)
ax.set_title('Effect of doubling sea salt emissions on surface temperature')

ax1=plt.subplot(122)
plt_violin(agregd,ax1)
ax1.set_title('Effect of doubling dust emissions on surface temperature')
ax1.set_ylabel('')

plt.show()
