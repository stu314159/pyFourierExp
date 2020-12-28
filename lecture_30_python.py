#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 10:21:01 2020

@author: Stu Blair

This is a Python implementation of an example problem from Lecture 30 of EM424.

The example is the solution of Laplace's Equation in Polar Coordinates.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import scipy.integrate as integrate

def ex1(theta):
    ret_val = 0.
    if (theta >= 0.) and (theta < np.pi/2.):
        ret_val = 1.;
    return ret_val;
           
            
def U_exp(theta,r,An,Bn,N):
    ret_val = 0.;
    for n in range(N):
        ret_val += (r**n)*(An[n]*np.cos(n*theta) + Bn[n]*np.sin(n*theta));
    return ret_val;
            
# parameters
c = 2.0; # radius 
N = 51; # number of Fourier Modes to be taken of the solution
# note that since Python uses zero-based indexing I will store Ao as A[0].
# B[0], by convention will be equal to zero. (yeah, that's a weird choice)
reltol = 1e-15;
# this example includes showing the students how to set custom error tolerances 

A = np.zeros(N,dtype=np.float64);
B = np.zeros_like(A); 

a_term, err = integrate.quad(lambda x: ex1(x),0.,2.*np.pi);
A[0] = a_term*(1./(2*np.pi));

for n in range(1,N):
    a_term, err = integrate.quad(lambda theta: ex1(theta)*np.cos(n*theta),
                                 0.,2.*np.pi,epsrel=reltol);
    A[n] = a_term*(1./((c**n)*np.pi));
    b_term, err = integrate.quad(lambda theta: ex1(theta)*np.sin(n*theta),
                                 0.,2.*np.pi,epsrel=reltol);
    B[n] = b_term*(1./((c**n)*np.pi));
    
# construct the solution in terms of the fourier coefficients
U = lambda theta,r: U_exp(theta,r,A,B,N);

NR = 100;
NT = 100;

R = np.linspace(0,c,NR);
T = np.linspace(0,2.*np.pi,NT);

RR, TT = np.meshgrid(R,T)

XX = np.multiply(RR,np.cos(TT));# is there a more simple syntax for this?
YY = np.multiply(RR,np.sin(TT));

UU = np.zeros_like(XX);
rows = XX.shape[0];
cols = XX.shape[1];

# if the below code is correct, it is a *lot* less convenient than the 
# MATLAB syntax.
for i in range(0,rows):
    for j in range(0,cols):
        UU[i,j] = U(TT[i,j],RR[i,j]);
        
fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(XX,YY,UU,cmap=cm.jet,
                       linewidth=0,antialiased=True)
ax.set_zlim(0,1.1)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.01f'));
ax.set_xlabel('X',fontsize=12,fontweight='bold');
ax.set_ylabel('Y',fontsize=12,fontweight='bold');
ax.set_zlabel('U',fontsize=12,fontweight='bold');
ax.view_init(15,235);

fig.colorbar(surf,shrink=0.5,aspect=10)
plt.title('Lecture 30 Example', fontsize=14, fontweight='bold')
plt.show()


