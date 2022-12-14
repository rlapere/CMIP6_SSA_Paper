########################
# Remy Lapere (03/11/2022)
# Figure 6 for Lapere et al., 2022 - CMIP6 SSaer
########################

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,BoundaryNorm
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy
import matplotlib as mpl
import math
import sys
from functions import *
mpl.rcParams.update({'font.size': 16})

u = 10
sst = 5

d1 = np.array([0.05,0.5,1,5])*2
d2 = np.array([0.05,0.5,1,5,10])*2
d3 = np.array([0.05,0.5,1,5,10,20])*2

ds = [d3,d2,d1]

def compute_sfs(diameters,u,tos):
    d_mean = []
    for k in range(len(diameters)-1):
        d_mean = np.append(d_mean,np.linspace(diameters[k],diameters[k+1],10))
    mm = 2200*4/3*math.pi*(d_mean/2*10**(-6))**3*10**(9)
    mass = (mm[1:]+mm[:-1])/2.0
    d_mean = d_mean[:-1]
    
    mo = monahan1986(d_mean,u)
    moflux = np.sum(mo*mass)

    sm = smithhar1998(d_mean,u)
    smflux = np.sum(sm*mass)

    go = gong2003(d_mean,u)
    goflux = np.sum(go*mass)
    
    ma = mahowald2006(d_mean,u)
    maflux = np.sum(ma*mass)*10**12

    mar = np.append(martensson2003(d_mean[d_mean<2.8],u,tos),monahan1986(d_mean[d_mean>=2.8],u))
    marflux = np.sum(mar*mass)

    sa = salter2015(d_mean,u,tos)
    saflux = np.sum(sa*mass)
    
    ja = jaegle2011(d_mean,u,tos)
    jaflux = np.sum(ja*mass)
    
    gr = grythe2014(d_mean,u,tos)
    grflux = np.sum(gr*mass)
    
    return [moflux,smflux,goflux,marflux,saflux,jaflux,maflux,grflux]

labs = ['MO86','SM98','GO03','MA03/\nMO86','SA15','JA11','MA06','GR14']
lb=['R$_{max}$=20'+u"\u03bc"+'m','R$_{max}$=10'+u"\u03bc"+'m','R$_{max}$=5'+u"\u03bc"+'m']

col=['teal','lightgray','steelblue']

ax0=plt.subplot(121)
i=0
for dd in ds:
    ous = compute_sfs(dd,u,sst)
    if i==0:
        maxous = ous
        minous = ous
    else:
        maxous = np.maximum(maxous,ous)
        minous = np.minimum(minous,ous)
    if i==1:
        plt.bar(range(len(ous)),ous,label=lb[i],color=col[i],edgecolor='k',linewidth=1,alpha=1,zorder=1)
    else:
        plt.plot(range(len(ous)),ous,markersize=16,linewidth=0,label=lb[i],color=col[i],markeredgecolor='k',marker='s',zorder=2)
    i=i+1

plt.vlines(range(len(ous)),maxous,minous,color='k')
plt.legend()
plt.xticks(range(len(ous)),labs,rotation=90)
plt.yscale('log')
plt.ylim(0.7,2000)
plt.ylabel('Mass flux ('+u"\u03bc"+'g$\,$m$^{-2}\,$s$^{-1}$)')
plt.title('a) Impact of cut-off size on mass emission flux')
plt.grid(axis='y',which='both')
ax0.set_axisbelow(True)

ax=plt.subplot(122)

dd=d2

uu=np.array([9,10,11])
lb=['U=9$\,$m/s','U=10$\,$m/s','U=11$\,$m/s']

i=0
for u_ in uu:
    ous = compute_sfs(dd,u_,sst)
    if i==0:
        maxous = ous
        minous = ous
    else:
        maxous = np.maximum(maxous,ous)
        minous = np.minimum(minous,ous)
    if i==1:
        plt.bar(range(len(ous)),ous,label=lb[i],color=col[i],edgecolor='k',linewidth=1)
    else:
        plt.plot(range(len(ous)),ous,markersize=16,linewidth=0,label=lb[i],color=col[i],markeredgecolor='k',marker='s')
    i=i+1

plt.vlines(range(len(ous)),maxous,minous,color='k')
plt.legend()
plt.xticks(range(len(ous)),labs,rotation=90)
plt.yscale('log')
plt.ylim(0.7,2000)
plt.title('b) Impact of wind speed on mass emission flux')
plt.grid(axis='y',which='both')
ax.set_axisbelow(True)

sss = [0,5,10]
colo=['seagreen','','turquoise']

i=0
for ss in sss:
    ous = compute_sfs(dd,10,ss)
    ous[0] = np.nan
    ous[1] = np.nan
    ous[2] = np.nan
    ous[6] = np.nan
    if i!=1:
        plt.plot(range(len(ous)),ous,markersize=20,linewidth=0,label=lb[i],color=colo[i],markeredgecolor='k',marker='*')
    i=i+1

plt.show()
