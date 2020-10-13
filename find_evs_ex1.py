#!/usr/bin/env python
"""
Written by: CAPT Stu Blair
Purpose: simple Fourier Series expansion exercises implmeneted in Python
"""

import numpy as np
import matplotlib.pyplot as plt


# define eigenfunctions
def y_n(alpha_n,x):
    return np.sin(alpha_n*x);
    
# define function with which to determin eigenvalues
def f(x):
    return x + np.tan(x);

# create a plot of f(x) to see what it looks like
a = 0.; b = 20.;
Nx= 1000;
X = np.linspace(a,b,Nx,endpoint=False)

plt.plot(X,f(X))
plt.axis([a,b,-20,20]);
plt.grid(True);
plt.xlabel('X',fontsize=16,fontweight='bold');
plt.ylabel('f(X)',fontsize=16,fontweight='bold');
plt.title('Non-linear function for eigenvalues',fontsize=16,fontweight='bold');
plt.show()