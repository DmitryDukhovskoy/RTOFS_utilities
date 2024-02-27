"""
	Find polynomial for an analytical function 
	mimicking observed depth-integrated transport
	Solve optimization problem: minimize error of polynomail
	interpolation given N mooring locations
"""
from IPython import get_ipython
get_ipython().magic('reset -sf')

import os
import numpy as np
from copy import copy
import importlib
import matplotlib.pyplot as plt
import sys

sys.path.append('/home/ddmitry/codes/MyPython/hycom_utils')
sys.path.append('/home/ddmitry/codes/MyPython/draw_map')
sys.path.append('/home/ddmitry/codes/MyPython')

from mod_utils_fig import bottom_text

def vflux(X):
	F = bet**alf*X**(alf-1.)*np.exp(-bet*X)
	return F

def lagr_polynom(Xp,Yp,xx):
	"""
	Lagrange polynomial
	estimate function at xx

	Pn(x) = sum(i to n): y(i)*Li(x)
	"""
	Np = Xp.shape[0]
	Pn = 0.
	for ii in range(Np):
		xi = Xp[ii]
		mi_x = 1.
		mi_xi = 1.
		for jj in range(Np):
			if jj == ii: 
				continue
			mi_xi = mi_xi*(xi-Xp[jj])  # mi(xi)
#			print('jj={0}, xi={1}, Xp={2}, mi_xi={3}'.format(jj,xi,Xp[jj],mi_xi))
			mi_x = mi_x*(xx-Xp[jj])  # mi(x)
#			print('xx={0}, mi_x={1}'.format(xx,mi_x))
	
		Li_x = mi_x/mi_xi
		Pn = Pn+Yp[ii]*Li_x

	return Pn

# Generate some continuous function 
# representing depth-integrated flux
X = np.arange(0,100.01,0.1)
alf = 3.
bet = 0.15
F = vflux(X) 


# Start with regularly distributed points:
N = 5
dx = (X[-1]-X[0])/N
Xp = np.arange(X[0],X[-1]+1.e-6,dx)
Yp = vflux(Xp)
Nord = N  # degree of polynomial

# Estimate at given points:
N = 100
dx = (X[-1]-X[0])/N
Xest = np.arange(X[0],X[-1]+1.e-6,dx)

Yest = []
for ip in range(N+1):
	yip = lagr_polynom(Xp,Yp,Xest[ip])
	Yest.append((yip))

# Exact solution:
Yex = vflux(Xest)
Err_inf = np.max(abs(Yex-Yest))/np.max(abs(Yex))



plt.ion()  # enables interactive mode
fig = plt.figure(1, figsize=(8,8))
fig.clf()
ax = plt.axes([0.1, 0.5,0.8,0.35])
ax.plot(X,F,label='analyt')
ax.plot(Xp,Yp,'r*',label='Interp Points')
ax.plot(Xest,Yest,label='Pn(x)')
ax.legend(loc='upper right')
plt.grid()

plt.title('Polynomial Pn(x), n={0}, Err_inf={1:.4f}'.format(Nord,Err_inf))
btx = 'polynom_obs.py'
bottom_text(btx)







