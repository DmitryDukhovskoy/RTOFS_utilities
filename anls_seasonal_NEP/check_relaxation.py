"""
  CHeck realxation schemes:
  Implicit scheme for the Newtonian term is 
  OPTIMAL NUDGING DATA ASSIMILATION SCHEME, Zou, Nvavon et al., 
Q. J. R. Meteorol. SOC. (1992), 118, pp. 116S1186 

  X(t+dt)-X_nn(t+dt) / 2*dt =  G(X_ref(t + dt) - X(t + dt)).

  X_nn - the tendency without the nudging term, X_nn(t + dt) being calculated first.

  In MOM6, the implicit scheme is written as:
          CS%var(m)%p(i,j,k) = I1pdamp * &
              (CS%var(m)%p(i,j,k) + CS%Ref_val(m)%p(k,c)*damp)

  where Ipdamp = 1/(1 + dt/tau), where tau - relaxation time, s

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import netCDF4
from netCDF4 import Dataset as ncFile
import importlib
import xarray
import yaml
from yaml import safe_load

PPTHN = []
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_read_hycom as mhycom
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_anls_seas as manseas


n=100
X = np.zeros((n))
X[0] = 1.2
dt = 200.
trlx = 7200. # relax. time, s
Irlx = 1./trlx  # s-1
I1pdamp = 1./(1. + dt*Irlx)
Xref = 2.85

for ii in range(1,n):
  X[ii] = I1pdamp*(X[ii-1] + dt*Irlx*Xref)


plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
ax1.plot(X)

