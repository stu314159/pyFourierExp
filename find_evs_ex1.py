#!/usr/bin/env python
"""
Written by: CAPT Stu Blair
Purpose: eigenvalue finding / Fourier Expansion process in Python
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize 
import scipy.integrate as integrate


# define eigenfunctions
def ef(alpha_n,x):
    return np.sin(alpha_n*x);

# construct Fourier expansion (this is crap)
def FF_exp(x,a_n,ev_n,N):    
    ret_val = 0.
    for i in range(N):
        ret_val += a_n[i]*ef(ev_n[i],x); # this really sucks
        
    return ret_val;
    
# define function with which to determine eigenvalues
def f(x):
    return x + np.tan(x);

# function that I want to expand in eigenfunctions
def f_exp(x):
    return np.exp(-x);

Nx= 1000;
# # create a plot of f(x) to see what it looks like
# a = 0.; b = 20.;
# X = np.linspace(a,b,Nx,endpoint=False)

# # plot the function to help develop algorithm for finding eigenvalues
# plt.plot(X,f(X))
# plt.axis([a,b,-20,20]);
# plt.grid(True);
# plt.xlabel('X',fontsize=16,fontweight='bold');
# plt.ylabel('f(X)',fontsize=16,fontweight='bold');
# plt.title('Non-linear function for eigenvalues',fontsize=16,fontweight='bold');
# plt.show()

N_evs = 35
roots = np.empty([N_evs,1]);

a = np.pi/2.+ 1e-6;
b = 3.*np.pi/2.-1e-6;
for i in range(N_evs):
    roots[i] = optimize.brentq(f,a,b);
    a = a + np.pi # reset lower bound
    b = b + np.pi # reset upper bound

# np.set_printoptions(precision=4);    
# print('The first %d roots:'%N_evs);
# print(roots);

# evs = roots**2.;
# print('The first %d eigenvalues:'%N_evs);
# print(evs);

#########################################################

# construct a Fourier Approximation using the eigenfunctions ef
a = 0; b = 1;
a_n = np.zeros_like(roots);

for i in range(N_evs):
    # get the next eigenfunction
    y_n = lambda x: ef(roots[i],x);
    
    # compute the Fourier Coefficient    
    num, err = integrate.quad(lambda x: f_exp(x)*y_n(x),a,b) 
    denom, err = integrate.quad(lambda x: y_n(x)**2,a,b)
    a_n[i] = num/denom;
        
   
FF = lambda x: FF_exp(x,a_n,roots,N_evs);


X = np.linspace(a,b,Nx)
plt.plot(X,f_exp(X),'-b',linewidth=2.0,label='f(x)')
plt.plot(X,FF(X),'--r',linewidth=2.0,label='FF(x)')
plt.grid(True);
plt.xlabel('X',fontsize=16,fontweight='bold');
plt.ylabel('f(X)',fontsize=16,fontweight='bold')
plt.title('Simple Fourier Expansion',fontsize=18,fontweight='bold')
plt.legend(loc='best');
plt.show();


    
