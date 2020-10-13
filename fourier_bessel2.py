#!/usr/bin/env python
"""
Written by: CAPT Stu Blair
Date: 12 Oct
Purpose: simple Fourier-Bessel series expansion to a function

This alternative implementation uses the more general jv function

"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.special import jn_zeros 
from scipy.special import jv 



def f_exp(x):
    return x;

def FB_exp(x,cn,ev,N):
    ret_val = 0.;
    for i in range(N):
        ret_val += cn[i]*jv(1,x*ev[i]);
    return ret_val;
    
# describe discrete domain
Nx = 1000;
a = 0.; b = 3.;
X = np.linspace(a,b,Nx);


# get eigenvalues
N = 35; # number of ev's.

# get N roots of J_1
k = jn_zeros(1,N);
# compute eigenvalues
ev = k/b;
# allocate an array for expansion coefficients
cn = np.zeros_like(ev);
for n in range(N):
    # compute all of the expansion coefficients
    num, err = integrate.quad(lambda x: f_exp(x)*jv(1,ev[n]*x)*x,a,b);
    denom, err = integrate.quad(lambda x: jv(1,ev[n]*x)*jv(1,ev[n]*x)*x,a,b);
    cn[n] = num/denom;
   
# make a convenient function capturing the FB expansion
FB = lambda x: FB_exp(x,cn,ev,N);

plt.plot(X,f_exp(X),'-b',linewidth=2.0,label='f(x)')
plt.plot(X,FB(X),'-.r',linewidth=2.0,label='FB(x)')
plt.grid(True);
plt.xlabel('X',fontsize=14,fontweight='bold');
plt.title('Fourier-Bessel Expansion',fontsize=16,fontweight='bold')
plt.legend()