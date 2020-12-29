#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 09:17:18 2020

@author: sblair

This is a Python implementation of an example problem from Lecture 31 of EM424.

The example is the solution of the Wave Equation in Polar Coordinates.

For this script I have implemented only the "ex2" initial conditions from the
MATLAB version of the example.
"""
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import scipy.integrate as integrate

from scipy.special import jn_zeros 
from scipy.special import j0 


# Parameters
N = 50; # number of modes
c = 1; # radius of the circular domain
a = 1; # "stiffness" parameter for the wave equation

# functions for initial conditions
def f(r):
    """
    Initial displacement for the wave equation.  The problem is assumed to 
    have radial symmetry so, for polar coordinates, the initial displacement
    is only a function of radial position.
    
    This is really a place-holder.  The example that I am implementing for
    lecture 31 has an initial displacement of 0.

    Parameters
    ----------
    r : float64 (or, whatever)
        radial position

    Returns
    -------
    float64 - initial displacement at a given radial position

    """
    return 0.

def g(r):
    """
    Initial velocity.  Radial symmetry so, for this problem, initial velocity
    is only a functin of radial position.

    Parameters
    ----------
    r : float64
        radial position.

    Returns
    -------
    float64 - initial velocity at a given radial position.

    """
    b = 0.2;
    v0 = 10.;
    retval = 0.;
    if r < b:
        retval = -v0;
    
    return retval;

# function to construct full solution
def U_exp(r,t,An,Bn,N):
    ret_val = 0.;
    for n in range(N):
        ret_val += j0(ev[n]*r)*(An[n]*np.cos(a*ev[n]*t) + 
                                Bn[n]*np.sin(a*ev[n]*t));
    return ret_val;


# allocate arrays
A = np.zeros(N,dtype=np.float64);
B = np.zeros_like(A); 

# get N roots of J_0
k = jn_zeros(0,N);
# compute eigenvalues
ev = k/c;

for n in range(N):
    ev_mag, err = integrate.quad(lambda r: r*j0(ev[n]*r)*j0(ev[n]*r),0,c);
    
    a_term, err = integrate.quad(lambda r: r*j0(ev[n]*r)*f(r),0,c);
    A[n] = a_term/ev_mag;
    
    b_term, err = integrate.quad(lambda r: r*j0(ev[n]*r)*g(r),0,c);
    b_term /= (a*ev[n]);
    B[n] = b_term/ev_mag;
    
# construct the solution in terms of the fourier coefficients
U = lambda r,t: U_exp(r,t,A,B,N);

# set-up geometric coordinates for plotting
NR = 20;
NT = 20;

R = np.linspace(0,c,NR);
T = np.linspace(0,2.*np.pi,NT);

RR, TT = np.meshgrid(R,T)

XX = np.multiply(RR,np.cos(TT));# is there a more simple syntax for this?
YY = np.multiply(RR,np.sin(TT));

UU = np.zeros_like(XX);
rows = XX.shape[0];
cols = XX.shape[1];

# since this is a time-dependent solution, we will want to make an animation.
# but first we will plot a solution at t=10 just to make sure everything
# is good-to-go.

t = 10.;

for i in range(0,rows):
    for j in range(0,cols):
        UU[i,j] = U(RR[i,j],t);
        
fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_wireframe(XX,YY,UU)

ax.set_zlim(-2.0,2.0)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.01f'));
ax.set_xlabel('X',fontsize=12,fontweight='bold');
ax.set_ylabel('Y',fontsize=12,fontweight='bold');
ax.set_zlabel('U',fontsize=12,fontweight='bold');
ax.view_init(20,235);
plt.title('Lecture 31 Example', fontsize=14, fontweight='bold')
plt.show()

# looks good, so I will try the time-dependent version
Tmax = 10
NTIME = 50
#T_time = np.linspace(0,Tmax,NTIME);

# below code based on example at: https://matplotlib.org/2.0.0/examples/mplot3d/wire3d_animation_demo.html

# function to generate "Z" values for the wireframe plot
def generate(RR,t,U):
    rows = RR.shape[0];
    cols = RR.shape[1];
    UU = np.zeros_like(RR);
    
    for i in range(0,rows):
        for j in range(0,cols):
            UU[i,j] = U(RR[i,j],t);
    return UU

# set-up the basic figure
fig = plt.figure();
ax = fig.add_subplot(111,projection='3d');
# use XX and YY from static plot
ax.set_zlim(-2.0,2.0)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.01f'));
ax.set_xlabel('X',fontsize=12,fontweight='bold');
ax.set_ylabel('Y',fontsize=12,fontweight='bold');
ax.set_zlabel('U',fontsize=12,fontweight='bold');
ax.view_init(20,235);
#plt.title('Lecture 31 Example', fontsize=14, fontweight='bold')

# begin plotting annimation
wframe = None

for t in np.linspace(0,Tmax,NTIME):
    if wframe:
        ax.collections.remove(wframe)
    UU = generate(RR,t,U);
    wframe = ax.plot_wireframe(XX,YY,UU)
    title_str = f"Lecture 31 example, t = {t:.2f}"
    plt.title(title_str,fontsize=14,fontweight='bold');
    plt.pause(0.01)
    
