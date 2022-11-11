###########
# Source function definitions for Figure 6
# Lapere et al., 2022 - CMIP6 SSaer
##########

import numpy as np
import math

def monahan1986(d,u):
    d=d/2
    B = (0.38-np.log10(d))/0.65
    W = 1.373*u**3.41
    mo = W*d**(-3)*(1.0+0.057*d**1.05)*10.0**(1.19*np.exp(-B**2))
    return mo


def smithhar1998(d,u):
    d=d/2
    A1 = 0.2*u**3.5
    A2 = 6.8*10.0**(-3)*u**3
    f1 = 1.5
    f2 = 1.0
    return A1*np.exp(-f1*np.log((d/3.0))**2)+A2*np.exp(-f2*np.log((d/30.0))**2)



def gong2003(d,u):
    d=d/2
    theta = 30
    A = 4.7*(1.0+theta*d)**(-0.017*d**(-1.44))
    B = (0.433-np.log10(d))/0.433
    W = 3.84*10**(-6)*u**(3.41)
    return W*3.6*10**5*d**(-A)*(1.0+0.057*d**(3.45))*10.0**(1.607*np.exp(-B**2))


def mahowald2006(d,u):
    S1 = 4.05*10**(-15)
    S2 = 4.52*10**(-14)
    S3 = 1.15*10**(-13)
    S4 = 1.20*10**(-13)
    W = u**3.41
    return S1*W*((d>=0.2) & (d<1.0))+S2*W*((d>=1.0) & (d<3.0))+S3*W*((d>=3.0) & (d<10.0))+S4*W*((d>=10.0) & (d<=40.0))


def martensson2003(d,u,t):
    W = 3.84*u**3.41*10**(-4)/100.0
    tw = t+273.15
    d = d*10**(-6)
    d_ = np.array([d**4,d**3,d**2,d,np.ones(len(d))])

    c0 = np.array([-2.576*10**35,5.932*10**28,-2.867*10**21,-3.003*10**13,-2.881*10**6])
    d0 = np.array([7.188*10**37,-1.616*10**31,6.791*10**23,1.829*10**16,7.609*10**8])
    A0 = d_[0]*c0[0]+d_[1]*c0[1]+d_[2]*c0[2]+d_[3]*c0[3]+d_[4]*c0[4]
    B0 = d_[0]*d0[0]+d_[1]*d0[1]+d_[2]*d0[2]+d_[3]*d0[3]+d_[4]*d0[4]
    phi0 = A0*tw+B0

    c1 = np.array([-2.452*10**33,2.404*10**27,-8.148*10**20,1.183*10**14,-6.743*10**6])
    d1 = np.array([7.368*10**35,-7.310*10**29,2.528*10**23,-3.787*10**16,2.279*10**9])
    A1 = d_[0]*c1[0]+d_[1]*c1[1]+d_[2]*c1[2]+d_[3]*c1[3]+d_[4]*c1[4]
    B1 = d_[0]*d1[0]+d_[1]*d1[1]+d_[2]*d1[2]+d_[3]*d1[3]+d_[4]*d1[4]
    phi1 = A1*tw+B1

    c2 = np.array([1.085*10**29,-9.841*10**23,3.132*10**18,-4.165*10**12,2.181*10**6])
    d2 = np.array([-2.859*10**31,2.601*10**26,-8.297*10**20,1.105*10**15,-5.8*10**8])
    A2 = d_[0]*c2[0]+d_[1]*c2[1]+d_[2]*c2[2]+d_[3]*c2[3]+d_[4]*c2[4]
    B2 = d_[0]*d2[0]+d_[1]*d2[1]+d_[2]*d2[2]+d_[3]*d2[3]+d_[4]*d2[4]
    phi2 = A2*tw+B2

    phi = phi0*((d>=0.02*10**(-6)) & (d<0.145*10**(-6)))+phi1*((d>=0.145*10**(-6)) & (d<0.419*10**(-6)))+phi2*((d>=0.419*10**(-6)) & (d<2.8*10**(-6)))
    return phi*W


def salter2015(d,u,t):
    A0,B0,C0,D0,d0,s0 = -5.2168*10**5,3.31725*10**7,-6.95275*10**8,1.0684*10**10,0.095,2.10
    A1,B1,C1,D1,d1,s1 = 0.0,7.374*10**5,-2.4803*10**7,7.7373*10**8,0.6,1.72
    A2,B2,C2,D2,d2,s2 = 0.0,1.4210*10**4,1.4662*10**7,1.7075*10**8,1.5,1.60

    N0 = A0*t**3+B0*t**2+C0*t+D0
    N1 = A1*t**3+B1*t**2+C1*t+D1
    N2 = A2*t**3+B2*t**2+C2*t+D2

    Fent = 3*10**(-8)*u**3.74

    flux0 = N0*Fent/(np.sqrt(2*math.pi)*np.log10(s0))*np.exp(-0.5*(np.log10(d)-np.log10(d0))**2/(np.log10(s0))**2)
    flux1 = N1*Fent/(np.sqrt(2*math.pi)*np.log10(s1))*np.exp(-0.5*(np.log10(d)-np.log10(d1))**2/(np.log10(s1))**2)
    flux2 = N1*Fent/(np.sqrt(2*math.pi)*np.log10(s2))*np.exp(-0.5*(np.log10(d)-np.log10(d2))**2/(np.log10(s2))**2)

    return flux0+flux1+flux2


def grythe2014(d,u,t):
    Tw = 0.3+0.1*t-0.0076*t**2+0.00021*t**3
    W = 235.0*u**3.5*np.exp(-0.55*(np.log(d/0.1))**2)+0.2*u**3.5*np.exp(-1.5*(np.log(d/3.0))**2)+6.8*u**3*np.exp(-1*(np.log(d/30.0))**2)
    return Tw*W

def jaegle2011(d,u,t):
    d=d/2
    theta = 30
    A = 4.7*(1.0+theta*d)**(-0.017*d**(-1.44))
    B = (0.433-np.log10(d))/0.433
    Tw = 0.3+0.1*t-0.0076*t**2+0.00021*t**3
    W = 3.84*10**(-6)*u**(3.41)
    fac=3.6*10**5*d**(-A)*(1.0+0.057*d**(3.45))*10.0**(1.607*np.exp(-B**2))
    return W*fac*Tw
                                               
def gongsalter(d,u,t):
    W = 3.84*10**(-6)*u**3.41
    d=d/2
    theta = 30
    A = 4.7*(1.0+theta*d)**(-0.017*d**(-1.44))
    B = (0.433-np.log10(d))/0.433
    fl= W*3.6*10**5*d**(-A)*(1.0+0.057*d**(3.45))*10.0**(1.607*np.exp(-B**2))
    fa = -8.75593*10.0**(-5)*t**3+5.56771*10**(-3)*t**2-0.11670*t+1.79321
    fc = 3.75294*10.0**(-2)*t+0.43706
    facc = fl*fa
    fcoa = fl*fc
    return np.append(facc[d<1],fcoa[d>=1])
