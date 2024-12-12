"""
  Investigate how accurate is to use climatologies for 
  removing seasonality
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import matplotlib.colors as colors
import random

PPTHN = '/home/Dmitry.Dukhovskoy/python'
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

random.seed

tau = 365.
A0  = 5.
mu      = 0.
sgm_a   = 0.2
sgm_eps = 0.6
sgm_tau = 5
ndays   = 365
nyrs = 20
t = np.arange(ndays)
y0 = 0
pdegr = 5 # polynomial degree

Yts = np.zeros((nyrs,ndays))
Eps_ts = np.zeros((nyrs,ndays))
for iyr in range(nyrs):
  a_rnd = random.uniform(0.8,1.2)
  t_rnd = random.normalvariate(mu, sgm_tau)
  eps   = np.random.normal(mu, sgm_eps, ndays)
  y = y0 + a_rnd*A0*np.sin(2*np.pi/(tau + t_rnd)*t) + eps
  Yts[iyr,:] = y
  Eps_ts[iyr,:] = eps

# Compute long-term mean:
Y_clim = np.nanmean(Yts, axis=0)

Eps_clim = np.zeros((nyrs,ndays))
Err_clim = np.zeros((nyrs))
for iyr in range(nyrs):
  eps_clim = Yts[iyr,:] - Y_clim
  Eps_clim[iyr,:] = eps_clim
  resid = eps_clim - Eps_ts[iyr,:]
  err = np.sqrt(np.sum(resid**2))
  Err_clim[iyr] = err

# Polynomial:
Eps_poly = np.zeros((nyrs,ndays))
Err_poly = np.zeros((nyrs))
for iyr in range(nyrs):
  y = Yts[iyr,:]
  Pcoef = np.polyfit(t, y, pdegr)
  Plnm  = np.poly1d(Pcoef)
  Pfit = Plnm(t)
  eps_poly = Yts[iyr,:] - Pfit
  Eps_poly[iyr,:] = eps_poly
  resid = eps_poly - Eps_ts[iyr,:]
  err = np.sqrt(np.sum(resid**2))
  Err_poly[iyr] = err

# 
# Running mean:
from mod_solver import runmn
# Filter N-point running mean
Eps_rmean = np.zeros((nyrs,ndays))
Err_rmean = np.zeros((nyrs))
tp1 = t.copy()
tp1 = np.append(tp1, tp1[-1]+1)  # running mean needs extra day 
for iyr in range(nyrs):
  y = Yts[iyr,:]
  Y_rmean = runmn(y, tp1, mnwnd=31)
  eps_rmean = Yts[iyr,:] - Y_rmean
  Eps_rmean[iyr,:] = eps_rmean
  resid = eps_rmean - Eps_ts[iyr,:]
  err = np.sqrt(np.sum(resid**2))
  Err_rmean[iyr] = err
  

def auto_corr(nlags, y1, y2):
  tlags = np.arange(-nlags,nlags+1)
  acorr = np.zeros((len(tlags)))*np.nan
  for ii in range(len(tlags)):
    ilag = tlags[ii]
    if ilag <= 0:
      ii1 = abs(ilag)
      y1s = y1[ii1:]
      na1 = len(y1s)
      y2s = y2[:na1]
    else:
      ii1 = ilag
      y2s = y2[ii1:]
      na1 = len(y2s)
      y1s = y1[:na1]

    sgm1 = np.std(y1s)
    sgm2 = np.std(y2s)
    corr = 1/ndays*np.sum((y1s-np.mean(y1s))*(y2s-np.mean(y2s)))/(sgm1*sgm2)

    acorr[ii] = corr
  
  return acorr, tlags

# Cross-correlation of not detrended time series:
# Random noise
y1 = Eps_ts[0,:]
y2 = Eps_ts[1,:]
# Noise with seasonal trend
y1tr = Yts[0,:]
y2tr = Yts[1,:]
nln = len(y1)
nlags = 100

acorr, tlags = auto_corr(nlags, y1, y2)
acorr_tr, _  = auto_corr(nlags, y1tr, y2tr)


# =====================
plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.24, 0.8, 0.7])
ax1.plot(Yts.transpose())
ax1.plot(Y_clim)
stl=f'y(t)=y0+eps_a*A*(cos(2*pi/(tau+eps_t)*t))+eps(~N(0,sgm))'
ax1.set_title(stl)
ax1.grid('on')

sinfo=f'eps_a~U(0.8, 1.2), A={A0:.2f}, tau={tau:.0f}, eps_t~N({mu:.1f}, {sgm_tau:.1f}), Nyears={nyrs}'
ax2 = plt.axes([0.1, 0.08, 0.8, 0.1])
ax2.text(0,0,sinfo)
ax2.axis('off')

btx = 'timeser_climat.py'
bottom_text(btx)

# Plot anomalies:
#iyr=10 use last year saved in memory
iyr=1
y = Yts[iyr,:]

Y_rmean = runmn(y, tp1, mnwnd=31)

Pcoef = np.polyfit(t, y, pdegr)
Plnm  = np.poly1d(Pcoef)
Pfit = Plnm(t)


# Show example for 1 year:
plt.clf()
ax1 = plt.axes([0.1, 0.24, 0.8, 0.7])
ax1.plot(y, color=[0, 0, 0])
ln1, = ax1.plot(Y_clim, label='Ensemble mean')
ln2, = ax1.plot(Pfit, label=f'P({pdegr})')
ln3, = ax1.plot(Y_rmean, label='running mean(30)')

stl=f'Year={iyr}'
ax1.set_title(stl)
ax1.grid('on')

sinfo=f'eps_a~U(0.8, 1.2), A={A0:.2f}, tau={tau:.0f}, eps_t~N({mu:.1f}, {sgm_tau:.1f}), Nyears={nyrs}'
ax2 = plt.axes([0.1, 0.08, 0.8, 0.1])
lgd = plt.legend(handles=[ln1,ln2,ln3], loc='upper right')
ax2.axis('off')

bottom_text(btx)


# Show anomalies for this year:
eps_ts = Eps_ts[iyr,:]
eps_clim = Eps_clim[iyr,:]
eps_poly = Eps_poly[iyr,:]
eps_rmean = Eps_rmean[iyr,:]


plt.clf()
ax1 = plt.axes([0.1, 0.24, 0.8, 0.7])
ln0, = ax1.plot(eps_ts, color=[0, 0, 0], label='True anomalies')
ln1, = ax1.plot(eps_clim, label='Ensemble mean')
ln2, = ax1.plot(eps_poly, label=f'P({pdegr})')
ln3, = ax1.plot(eps_rmean, label='running mean(30)')

stl=f'Year={iyr}, recovered anomalies'
ax1.set_title(stl)
ax1.grid('on')

ax2 = plt.axes([0.1, 0.08, 0.8, 0.1])
lgd = plt.legend(handles=[ln0,ln1,ln2,ln3], loc='upper right')
ax2.axis('off')

bottom_text(btx)

# Show error:
plt.clf()

ax1 = plt.axes([0.1, 0.24, 0.8, 0.7])
ln1, = ax1.plot(Err_clim, label='Ensemble mean')
ln2, = ax1.plot(Err_poly, label=f'P({pdegr})')
ln3, = ax1.plot(Err_rmean, label='running mean(30)')

stl=f'sqrt(||res||2) for recovered anomalies'
ax1.set_title(stl)
ax1.grid('on')
ax1.set_xlabel('Years')

ax2 = plt.axes([0.1, 0.08, 0.8, 0.1])
lgd = plt.legend(handles=[ln1,ln2,ln3], loc='upper right')
ax2.axis('off')

bottom_text(btx)


# Plot correlation:
plt.clf()

ax1 = plt.axes([0.1, 0.24, 0.8, 0.7])
ln1, = ax1.plot(tlags, acorr, label='Random noise')
ln2, = ax1.plot(tlags, acorr_tr, label='Random noise + seasonal')


stl='Cross-correlation for random time series with/without seasonal signal'
ax1.set_title(stl)
ax1.grid('on')
ax1.set_xlabel('Lags')
ax1.set_ylim(-0.2, 1.)


ax2 = plt.axes([0.1, 0.08, 0.8, 0.1])
lgd = plt.legend(handles=[ln1,ln2], loc='upper right')
ax2.axis('off')

bottom_text(btx)



