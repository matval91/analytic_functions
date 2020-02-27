#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to compute relaxation rate against feeding rate for NNB on JT60SA
I am using formulas at page 31 of NRL formulary
Created on Wed Feb 12 09:54:41 2020

@author: vallar
"""
import numpy as np
import matplotlib.pyplot as plt
import sympy
import scipy.constants as const
## Defining constants
qe = const.e; mp=const.proton_mass;
A=2; Z=1 # A and Z of beam ions
me = const.electron_mass
lnlambda = 17
kb = const.k #in m2 kg s-2 K-1
epsilon0 = const.epsilon_0 #in A2⋅s4⋅kg−1⋅m−3
## Defining variables
Eev = 500e3 # Energy in keV
Ej = Eev*qe # energy in Joule
va = np.sqrt(2*Ej/(A*mp)) # velocity of incoming fast ions

I_NB = 10e6/Eev #Power over DV!
rate_NB = I_NB/(Z*qe) #number of ions per second

#
#Eev=85
#Ej = Eev*qe # energy in Joule
#va = np.sqrt(2*Ej/(A*mp)) # velocity of incoming fast ions
#I_NB = 24e6/Eev #Power over DV!
#rate_NB = I_NB/(Z*qe) #number of ions per second

v=np.linspace(va*0.2, va*1.1, 20)
E = v*v*A*mp*0.5
E = E/qe
nue = np.zeros(np.size(v))
nui = np.zeros(np.size(v))
target='e'
## now considering targetting on electrons
n_beta = 8.7e19 #line-averaged density of scenario 3 from J.Garcia's paper
Te_ev = 4.1e3 #line-averaged te of scenario 3 from J.Garcia's paper
Te_j = Te_ev*qe
ve_t = np.sqrt(2*Te_j/(me))#sqrt(2kT/m)
gamma = (Z*qe)**2*(qe)**2*lnlambda/(8*np.pi*epsilon0**2)
nu_0 = (1+A*mp/me)*4*gamma*n_beta/((A*mp)**2*ve_t**3)

#######
#OTHER FORMULA OF COLLISION FREQUENCY
#nu = (4*pi*ne*Zi*e**4*lnZ)/((4*pi*eps0)**2me**2*V**3)
nu=np.zeros(np.size(v))

#for i,vel in enumerate(v):
nu = (4*const.pi*n_beta*qe**4*lnlambda)/((4*const.pi*epsilon0)**2*me**2*v**3)
#
#    
#    print('doing iteration {:d} over {:d}'.format(i, np.size(v)))
#    x=vel/ve_t
#    eta = sympy.Symbol('eta')
#    phi = sympy.erf(eta)
#    dphi = sympy.diff(phi,eta)
#    phi_x = float(phi.evalf(subs={eta:x}))
#    dphi_x = float(dphi.evalf(subs={eta:x}))
#    psi = 1/(2*x**2)*(phi_x-x*dphi_x)
#    nue[i] = nu_0*psi/x
## target on ions
##n_beta = 7e-19 #line-averaged density of scenario 3 from J.Garcia's paper
#Ti_ev = 4.2e3 #line-averaged ti of scenario 3 from J.Garcia's paper
#Ti_j = Ti_ev*qe
#vi_t = np.sqrt(2*Ti_j/(A*mp))#sqrt(2kT/m)
#gamma = (Z*qe)**2*(qe)**2*lnlambda/(8*np.pi*epsilon0**2)
#nu_0 = (1+1)*4*gamma*n_beta/((A*mp)**2*vi_t**3)
#
#for i,vel in enumerate(v):
#    print('doing iteration {:d} over {:d}'.format(i, np.size(v)))
#    x=vel/ve_t
#    # t=sympy.Symbol('t')
#    eta = sympy.Symbol('eta')
#    phi = sympy.erf(eta)
#    dphi = sympy.diff(phi,eta)
#    phi_x = float(phi.evalf(subs={eta:x}))
#    dphi_x = float(dphi.evalf(subs={eta:x}))
#    # #psi = 2/np.sqrt(np.pi)*sympy.integrate(sympy.exp(-t)*t**(0.5), (t, 0, x))
#    # #psi = float(psi.evalf())
#    psi = 1/(2*x**2)*(phi_x-x*dphi_x)
#    nui[i] = nu_0*psi/x



plot=True
if plot:
    f=plt.figure()
    ax=f.add_subplot(111)
#    ax.plot(E*1e-3, nue, lw=2.3, color='b', label='e')
#    ax.plot(E*1e-3, nui, lw=2.3, color='r', label='i')
#    ax.plot(E*1e-3, nui+nue, lw=2.3, color='k', label='tot')
    ax.plot(E*1e-3, nu, lw=2.3, color='k', label='tot')

    ax.axvline(Eev*1e-3)
    ax.plot([min(E)*1e-3, max(E)*1e-3], [rate_NB, rate_NB], 'k--', label='NB')
    ax.legend()
    ax.set_xlabel(r'E [keV]')
    ax.set_ylabel(r'$\nu$ [1/s]')
#    
#    f2=plt.figure(); ax2=f2.add_subplot(111)
#    ax2.plot(E*1e-3, -(nui+nue)*v)
#    ax2.set_xlabel(r'E')
#    ax2.set_ylabel(r'dv/dt')