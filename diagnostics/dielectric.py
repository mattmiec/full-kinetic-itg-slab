# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 15:33:49 2018

@author: mtmiec
"""

from math import pi,sqrt
from cmath import exp
from scipy.special import wofz,ive
#Exponentially scaled modified Bessel function of the first kind
#ive(v, z) = iv(v, z) * exp(-abs(z.real))

def plasmaDisp(z):
    """
    returns tuple (Z,Zp,Zpp) of plasma dispersion function Z,Z',Z''
    plasma dispersion functoin is defined as
    i*sqrt(pi)*exp(-z^2)*erfc(-i*z)
    """
    Z = 1j*sqrt(pi)*wofz(z)
    Zp = -2*(1+z*Z)
    Zpp = -2*(Z+z*Zp)
    return (Z, Zp, Zpp)

def dielectricFunc(omega,theta,kx,ky,kpar,kapn,kapt):
    """
    returns dielectric function for given parameters
    omega is complex frequency normalized to omega_i, expressed as 2 element array of real, imag parts
    theta is ratio q_i T_e / e T_i
    kx,ky,kpar are wavenumbers normalized to rho_i
    kapn,kapt are inverse gradient scale lengths normalized to rho_i
    """
    
    #convert to complex
    omega = (omega[0]+1j*omega[1])/sqrt(2)/kpar
    
    b = kx**2 + ky**2
    (plz,plzp,plzpp) = plasmaDisp(omega)
    gm0 = ive(0,b)
    gm1 = ive(1,b)

    tol = 1e-12
    
    tgrad = -theta*ky*kapt*(b*(gm0-gm1)*plz-.25*gm0*plzpp)/kpar/sqrt(2)
    ngrad = theta*ky*kapn*gm0*plz/kpar/sqrt(2)
    iaw = -.5*theta*gm0*plzp
    pol = 0.
    
    n = 0 
    gmnm1 = gm0
    gmn = gm1
    
    flag = 0
    nmax = 10000
    
    while flag == 0:
        n=n+1
        gmnp1 = ive(n+1,b)
        zetp = omega+n/sqrt(2)/kpar
        zetm = omega-n/sqrt(2)/kpar

        (pplz,pplzp,pplzpp) = plasmaDisp(zetp)
        (mplz,mplzp,mplzpp) = plasmaDisp(zetm)
    
        tgradc = -theta*ky*kapt*(.5*b*(2*gmn-gmnp1-gmnm1)*(pplz+mplz)-.25*gmn*(pplzpp+mplzpp))/kpar/sqrt(2)
        ngradc = theta*ky*kapn*gmn*(pplz+mplz)/kpar/sqrt(2)
        iawc = -.5*theta*gmn*(pplzp+mplzp)
        polc = -theta*gmn*n*(pplz-mplz)/kpar/sqrt(2)

        tgrad = tgrad + tgradc
        ngrad = ngrad + ngradc
        iaw = iaw + iawc
        pol = pol + polc
        
        gmnm1 = gmn
        gmn = gmnp1
        
        if abs(tgradc)/abs(tgrad) > tol:
            flag = 0
        elif abs(ngradc)/abs(ngrad) > tol:
            flag = 0
        elif abs(iawc)/abs(iaw) > tol:
            flag=0
        elif abs(polc)/abs(pol) > tol:
            flag=0
        else:
            flag=1

        if n > nmax:
            flag=1
            print('Warning: exceeded max iteratations')

    ep = 1 + tgrad + ngrad + iaw + pol
    return [ep.real,ep.imag]
    
def dielectricFunc_kpar0_ky0(omega,theta,kx,alpha):
    """
    returns dielectric function (* alpha kx^2) for given parameters
    neglecting gradients, with k_par = 0 and k_y = 0
    omega is complex frequency normalized to omega_i, expressed as 2 element array of real, imag parts
    theta is ratio q_i T_e / e T_i
    kx is wavenumber normalized to rho_i
    """
    
    #convert to complex
    omega = (omega[0]+1j*omega[1])
    b = kx**2
    tol = 1e-8    
    pol = 0.    
    n = 0
    gm0 = ive(0,b)
    gm1 = ive(1,b)
    gmnm1 = gm0
    gmn = gm1    
    flag = 0
    nmax = 10000
    
    while flag == 0:
        n=n+1
        gmnp1 = ive(n+1,b)
        zetp = omega+n
        zetm = omega-n

        polc = theta*gmn*n*(1/zetp-1/zetm)

        pol = pol + polc
        
        gmnm1 = gmn
        gmn = gmnp1
        
        if abs(polc)/abs(pol) > tol:
            flag=0
        else:
            flag=1

        if n > nmax:
            flag=1
            print('Warning: exceeded max iteratations')

    ep = alpha*b + pol
    return [ep.real,ep.imag]


def DKEdielectricFunc(omega,theta,kx,ky,kpar,kapn,kapt,alpha,memi):
    """
    returns dielectric function for given parameters
    omega is complex frequency normalized to omega_i, expressed as 2 element array of real, imag parts
    theta is ratio q_i T_e / e T_i
    kx,ky,kpar are wavenumbers normalized to rho_i
    kapn,kapt are inverse gradient scale lengths normalized to rho_i
    """
    
    #convert to complex
    omega = (omega[0]+1j*omega[1])/sqrt(2)/kpar

    #electron part
    omega_e = omega*(memi**.5)
    (plz,plzp,plzpp) = plasmaDisp(omega_e)
    electro = -.5*plzp    
    
    b = kx**2 + ky**2
    (plz,plzp,plzpp) = plasmaDisp(omega)
    gm0 = ive(0,b)
    gm1 = ive(1,b)

    tol = 1e-8
    
    tgrad = -theta*ky*kapt*(b*(gm0-gm1)*plz-.25*gm0*plzpp)/kpar/sqrt(2)
    ngrad = theta*ky*kapn*gm0*plz/kpar/sqrt(2)
    iaw = -.5*theta*gm0*plzp
    pol = 0.
    
    n = 0 
    gmnm1 = gm0
    gmn = gm1
    
    flag = 0
    nmax = 10000
    
    while flag == 0:
        n=n+1
        gmnp1 = ive(n+1,b)
        zetp = omega+n/sqrt(2)/kpar
        zetm = omega-n/sqrt(2)/kpar

        (pplz,pplzp,pplzpp) = plasmaDisp(zetp)
        (mplz,mplzp,mplzpp) = plasmaDisp(zetm)
    
        tgradc = -theta*ky*kapt*(.5*b*(2*gmn-gmnp1-gmnm1)*(pplz+mplz)-.25*gmn*(pplzpp+mplzpp))/kpar/sqrt(2)
        ngradc = theta*ky*kapn*gmn*(pplz+mplz)/kpar/sqrt(2)
        iawc = -.5*theta*gmn*(pplzp+mplzp)
        polc = -theta*gmn*n*(pplz-mplz)/kpar/sqrt(2)

        tgrad = tgrad + tgradc
        ngrad = ngrad + ngradc
        iaw = iaw + iawc
        pol = pol + polc
        
        gmnm1 = gmn
        gmn = gmnp1
        
        if abs(tgradc)/abs(tgrad) > tol:
            flag = 0
        elif abs(ngradc)/abs(ngrad) > tol:
            flag = 0
        elif abs(iawc)/abs(iaw) > tol:
            flag=0
        elif abs(polc)/abs(pol) > tol:
            flag=0
        else:
            flag=1

        if n > nmax:
            flag=1
            print('Warning: exceeded max iteratations')

    ep = alpha*b + tgrad + ngrad + iaw + pol + electro
    return [ep.real,ep.imag]
