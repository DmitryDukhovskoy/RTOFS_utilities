"""
  Plot ice thickness NSIDC data
  satellites
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
from copy import copy
import matplotlib.colors as colors
from yaml import safe_load

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
import mod_plot_xsections as mxsct
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
#import mod_valid_utils as mvutil
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_anls_seas as manseas
import mod_utils_ob as mutob
importlib.reload(mutob)

pthdata = '/work/Dmitry.Dukhovskoy/data/NSIDC_icethckn'
flithkn = 'ers1_seaice_thickness_mean_march_1993to2001.nc'
dflithkn = os.path.join(pthdata, flithkn)

dset = xarray.open_dataset(dflithkn)
LAT  = dset['lat'].data
LON  = dset['lon'].data
HI   = dset['thickness'].data

clrmp = mclrmps.colormap_ice_thkn()
clrmp.set_bad(color=[0.2, 0.2, 0.2])
rmin = 0.
rmax = 5.

Nclrs = clrmp.N
cval  = np.linspace(rmin, rmax, Nclrs)

sttl = f"NSIDC {flithkn}" 
# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
#m = Basemap(width=3300*1.e3,height=3300*1.e3, resolution='l',\
#            projection='stere', lat_ts=55, lat_0=62, lon_0=-175)

m = Basemap(projection='npstere',boundinglat=50,lon_0=0,resolution='l')

LON = np.where(LON<-900, np.nan, LON)
LAT = np.where(LAT<-900, np.nan, LAT)

xR, yR = m(LON, LAT)
PMsk = ( (xR > 1e20) | (yR > 1e20) )
#HI[PMsk] = 0.
#xR[PMsk] = 1.e30
#yR[PMsk] = 1.e30

VMsk = ( (xR < 1.e20) & (yR < 1.e20) )
J, I = np.where( (~np.isnan(HI))  &  (HI>0.) )

#img = m.pcolormesh(xR, yR, HI, cmap=clrmp, vmin=rmin, vmax=rmax)
vclrs = np.zeros((len(J),3))
for ii in range(len(J)):
  j0 = J[ii]
  i0 = I[ii]
  hh = HI[j0,i0]
  if hh > cval[-1]:
    vclrs[ii,:] = clrmp(Nclrs)[:3]
  else:
    iclr = max(np.where(cval <= hh)[0])
    vclrs[ii,:] = clrmp(iclr)[:3] 

msize = np.zeros((len(J))) + 0.1

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

m.scatter(xR[J,I], yR[J,I], c=vclrs, s=8)



