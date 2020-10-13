#!/usr/bin/env python
"""
Written by: CAPT Stu Blair
Purpose: eigenvalue finding / Fourier Expansion process in Python
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize 


# define eigenfunctions
def y_n(alpha_n,x):
    return np.sin(alpha_n*x);
    
# define function with which to determine eigenvalues
def f(x):
    return x + np.tan(x);

# function that I want to expand in eigenfunctions
def f_exp(x):
    return np.exp(-x);

# create a plot of f(x) to see what it looks like
a = 0.; b = 20.;
Nx= 1000;
X = np.linspace(a,b,Nx,endpoint=False)

# plot the function to help develop algorithm for finding eigenvalues
plt.plot(X,f(X))
plt.axis([a,b,-20,20]);
plt.grid(True);
plt.xlabel('X',fontsize=16,fontweight='bold');
plt.ylabel('f(X)',fontsize=16,fontweight='bold');
plt.title('Non-linear function for eigenvalues',fontsize=16,fontweight='bold');
plt.show()

N_evs = 8
roots = np.empty([N_evs,1]);

a = np.pi/2.+ 1e-6;
b = 3.*np.pi/2.-1e-6;
for i in range(N_evs):
    roots[i] = optimize.brentq(f,a,b);
    a = a + np.pi # reset lower bound
    b = b + np.pi # reset upper bound

np.set_printoptions(precision=4);    
print('The first %d roots:'%N_evs);
print(roots);

evs = roots**2.;
print('The first %d eigenvalues:'%N_evs);
print(evs);
    
