"""
	Find polynomial for observed depth-integrated transport
	picewise polynomial: try linear, quadr, cubic
	then Hermite polynomial
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

plt.close('all')
rVFlx = False   # rVFlx - read/plot Vol Flux, otherwise - FW flux

from mod_utils_fig import bottom_text
plt.ion()  # enables interactive mode


def read_vflux(yr,strnm):
	pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_straits/'
	flxfile = pthout+'Fluxes_'+strnm+'_{0}.dat'.format(yr)
	print('Reading ',flxfile)

	fid = open(flxfile,'rb')
	fid.seek(0)
	dmm = np.fromfile(fid,dtype='>i',count=1)
	nflds = dmm[0]
	dmm = np.fromfile(fid,dtype='>i',count=1)
	npnts = dmm[0]
	dL = np.fromfile(fid,dtype='>f',count=npnts)
	Vflx = np.fromfile(fid,dtype='>f',count=npnts)
	Tflx = np.fromfile(fid,dtype='>f',count=npnts)
	FWflx = np.fromfile(fid,dtype='>f',count=npnts)
	fid.close()

	Vflx = 1.e-6*Vflx  # Sv
	return dL, Vflx, Tflx, FWflx

def lagr_polynom(Xp,Yp,xx):
	"""
	Lagrange polynomial
	estimate function at xx

	Pn(x) = sum(i to n): y(i)*Li(x)

	Note that Lagrange polynom. parameters (Li_x) can be used for any
	Y, so no need to recompute parameters as long as grid points are the same

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

def plot_flx(X,F,Xp,Yp,Xest,Peqd,ERR,n1,n2,stl=' title ',xlbl='East Dist, km', ylbl='Flux, Sv'):
	plt.ion()  # enables interactive mode
	fig = plt.figure(figsize=(8,8))
	fig.clf()
	ax = plt.axes([0.1, 0.6,0.85,0.35])
	ax.plot(X,F,label='Obs')
	ax.plot(Xp,Yp,'r*',label='Interp Points')
	ax.plot(Xest,Peqd,label='Pn(x)')
	ax.legend(bbox_to_anchor=(0.98,-0.1))
	ax.set_xlabel(xlbl)
	ax.set_ylabel(ylbl)
	plt.grid()
	ax.set_title(stl)

# Plot Error
	xr = np.arange(n1,n2+1)
	ax2 = plt.axes([0.1,0.1,0.8,0.35])
	ax2.plot(xr,ERR,'.-')
	ax2.set_yscale("log")
	ax2.set_xlabel('# of intrp points')
	ax2.set_ylabel('||Err||_{inf}')
	ax2.set_title('Error vs N. of interpolation points')
	plt.grid()

	btx = 'polynom_obs.py'
	bottom_text(btx)



yr = 2019
strnm = 'FramStr'

dL, Vflx, Tflx, FWflx = read_vflux(yr,strnm)

if rVFlx:
	Flx = Vflx
	FlxNm = 'VolFlux'
else:
	Flx = FWflx
	FlxNm = 'FWFlux'

# Construct distance array in 100 km
XX = np.cumsum(dL*1.e-3) 
XX = XX-dL[0]*1.e-3

# Equidistant points:
Xeqd = XX


# Start with regularly distributed points:
nP = XX.shape[0]
ERR_eqd = []
ERR_chb = []
n1 = 4
n2 = 26
iplt = 9
if iplt>n2:
	iplt = n2

for N in range(n1,n2+1):
	dx = (nP-1)/(N-1)
	iXp = np.round(np.arange(0,nP-0.1,dx)).astype(int)
	iXp = list(iXp)

	Yp = Flx[iXp]
	Xp = XX[iXp]

# Chebyshev 1st kind  points:
	kk = np.arange(N+1)
	Ichb = np.cos(np.pi*(2.*kk+1.)/(2.*N+2.))  # on -1 < x < +1
# Convert to dist:
	Ichb = (1.+Ichb)/2.0
	Ichb = np.sort(Ichb)
# Find closest indices: - may not be very accurate Chebyshev points:
	iXch = np.round(Ichb*(nP-1)).astype(int)
	Ychb = Flx[iXch]
	Xchb = XX[iXch]	


	Peqd = []  # estimated flux for equidistant polynomial points
	Pchb = []  # estimated flux for Chbyshev pnts
	for ip in range(nP):
		yip = lagr_polynom(Xp,Yp,XX[ip])
		Peqd.append((yip))

		yip = lagr_polynom(Xchb,Ychb,XX[ip])
		Pchb.append((yip))

	if N == iplt:
		Peqd_plt = Peqd.copy()
		Xeqd_plt   = Xp.copy()
		Yeqd_plt   = Yp.copy()

		Pchb_plt = Pchb.copy()
		Xchb_plt = Xchb.copy()
		Ychb_plt = Ychb.copy()
		

# Exact solution:
	Yex = Flx
	Err_inf1 = np.max(abs(Yex-Peqd))/np.max(abs(Yex))
	ERR_eqd.append((Err_inf1))

	Err_inf2 = np.max(abs(Yex-Pchb))/np.max(abs(Yex))
	ERR_chb.append((Err_inf2))
	print('N pnts={0}, Err_eqd={1:6.3f} Err_chb={2:6.3f}'.format(N,Err_inf1,Err_inf2))

err = ERR_eqd[iplt]
stl = 'Fram {2}, HYCOM0.08, Eqdist pnts, Polynomial Pn(x), n={0}, Err_inf={1:.4f}'.\
			format(iplt,err,FlxNm)
plot_flx(XX,Flx,Xeqd_plt,Yeqd_plt,XX,Peqd_plt,ERR_eqd,n1,n2,stl,ylbl='VFlux, Sv')

err = ERR_chb[iplt]
stl = 'Fram {2}, HYCOM0.08, Cheb.pnts, Polynomial Pn(x), n={0}, Err_inf={1:.4f}'.\
			format(iplt,err,FlxNm)
plot_flx(XX,Flx,Xchb_plt,Ychb_plt,XX,Pchb_plt,ERR_chb,n1,n2,stl,ylbl='VFlux, Sv')

# Plot last polynom for Chebyshev
err = ERR_chb[-1]
stl = 'Fram {2}, HYCOM0.08, Cheb.pnts, Polynomial Pn(x), n={0}, Err_inf={1:.4f}'.\
			format(n2,err,FlxNm)
plot_flx(XX,Flx,Xchb,Ychb,XX,Pchb,ERR_chb,n1,n2,stl,ylbl='VFlux, Sv')








