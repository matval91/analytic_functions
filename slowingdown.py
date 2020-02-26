#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 16:12:55 2020

@author: vallar
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy
## Defining constants
qe = 1.602e-19; mp=1.66e-27;
A=2; Z=1 # A and Z of beam ions
me = 9.10938356e-31
lnlambda = 17
kb = 1.38064852e-23 #in m2 kg s-2 K-1
epsilon0 = 8.85e-12 #in A2⋅s4⋅kg−1⋅m−3
## Defining variables
Eev = 500e3 # Energy in keV
Ej = Eev*qe # energy in Joule
va = np.sqrt(2*Ej/(A*mp)) # velocity of incoming fast ions

v=np.linspace(va*0.2, va*1.1, 20)
E = v*v*A*mp*0.5
E = E/qe
## now considering targetting on electrons
n_beta = 8.7e19 #line-averaged density of scenario 3 from J.Garcia's paper
Te_ev = 4.1e3 #line-averaged te of scenario 3 from J.Garcia's paper
Te_j = Te_ev*qe
ve_t = np.sqrt(2*Te_j/me)#sqrt(2kT/m)

EC=100e3;
E=sympy.Symbol('E')
aaa=sympy.sqrt(E)/((qe*EC)**sympy.Rational(3,2)+E**sympy.Rational(3,2))
energies = np.linspace(0e3, 520e3, 200)
energies *=qe
sd=np.zeros(np.size(energies))
for i,el in enumerate(energies):
    try:
        t = sympy.integrate(aaa, (E, energies[i-1], el))
    except:
        sd[i]=0;
        continue
    
    sd[i]=float(t.evalf())
factor=6.27e8*A*Te_ev**1.5/(Z**2*n_beta*lnlambda)


f=plt.figure(); ax=f.add_subplot(111)
ax.plot(energies/qe, sd, 'k+', lw=2.3)
ax.axvline(85e3, color='r')
ax.axvline(500e3, color='g')
ax.set_xlabel(r'E [J]')
ax.set_ylabel(r'SD time (AU)')
ax.set_ylim([0, 1])
ax.grid('on')
f.tight_layout()