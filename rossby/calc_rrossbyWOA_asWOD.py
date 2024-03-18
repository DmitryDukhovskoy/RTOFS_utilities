# The code is similar to calc_rrossby_WOD.py
# but computes Rossby radius for 1 profile from the WOA
# specified by lon/lat (or grid indices) 
#
# For testing and comparison with WOA
# 
# 1st baroclinic Rossby radius 
#
# Use bottom topography for depths deeper than WOA last depth level (-5500m)
#  Topography interpolated from ETOPO2
#
# Calculate N2 following Chelton 1996
# N2 is calculated at the middle of the depth interval
# density is adiabatically adjusted to the mid-grid depth
#
import os
import numpy as np
import matplotlib.pyplot as plt
#import torch
import sys
import pdb
import netCDF4
import importlib
from netCDF4 import Dataset as ncFile
import timeit

sys.path.append('/data/home004/dmitry.dukhovskoy/python_codes/MyPython')
sys.path.append('/data/home004/dmitry.dukhovskoy/python_codes/read_wod')

import mod_utils_fig
# importlib.reload(mod_utils_fig)
from mod_utils_fig import bottom_text
from mod_utils_fig import plot_fld2D

# ========================================
# Specify location and depths where T/S data are "observed"
lon_plt = 137.
lat_plt = 0.0

Zintrp = np.array([-1.5,  -25. ,  -50. ,  -75. , -100. , -125. , -150. , -200. ,
       -250. , -300. , -500. , -750.])
#Zobs = []  # keep empty if no subsetting is needed
# ========================================

grd=0.25
if grd==0.25:
  cgrd=4
woa='woa18'
seas=13    # season: 1-12 - monthly, 13-16 - spring, summer, ...

pthout = '/data/ncei2/work/dmitry.dukhovskoy/data_output/'

btx = 'calc_rrossbyWOA_asWOD.py'
urlT='https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/temperature/decav/{0}/'.\
      format(grd)
urlS='https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/salinity/decav/{0}/'.\
      format(grd)

tfnm='{0}_decav_t{1:02d}_{2:02d}.nc'.format(woa,seas,cgrd)
sfnm='{0}_decav_s{1:02d}_{2:02d}.nc'.format(woa,seas,cgrd)

#

def read_field(furl,varnm):
  print("Reading {1} from {0}".format(furl,varnm))
  nc=ncFile(furl)
# lookup a variable
  dmm0 = nc.variables[varnm][:].squeeze()
  dmm = np.copy(dmm0)
  return dmm

def read_1prof(furl,varnm,j0,i0):
  """
    Read 1 vertical profile for specified location
    Note indices in netcdf run 0,... - same as in Python
  """
  nc=ncFile(furl)
  dmm0 = nc.variables[varnm][0,:,j0,i0].squeeze()
  dmm = np.copy(dmm0)
  return dmm

def get_dim(furl):
  nc=ncFile(furl)
  kdm = nc.dimensions['depth'].size
  jdm = nc.dimensions['lat'].size
  idm = nc.dimensions['lon'].size

  return idm, jdm, kdm


iz=0 
tvar='t_an'
svar='s_an'

furl=urlT+tfnm
furlT = furl
#T=read_field(furl,tvar)
#T[T>1.e10]=np.nan
#kdm=T.shape[0]
#jdm=T.shape[1]
#idm=T.shape[2]

idm, jdm, kdm = get_dim(furl)

ZZ=read_field(furl,'depth')
ZZ=-abs(ZZ)
lat=read_field(furl,'lat')
lon=read_field(furl,'lon')

LON=np.zeros((jdm,idm))
LAT=np.zeros((jdm,idm))
for ii in range(idm):
  LAT[:,ii]=lat

for jj in range(jdm):
  LON[jj,:]=lon


furlS=urlS+sfnm
#S=read_field(furlS,svar)
#S[S>1.e10]=np.nan

# Land-sea mask
#aa=T[iz,:,:].squeeze()
#ss=S[iz,:,:].squeeze()
#LMsk=aa.copy();
#LMsk[np.isfinite(aa)]=1.
#LMsk[~np.isfinite(aa)]=0.

# ------------
# Read TOPO
# -------------
pthout = '/data/ncei2/work/dmitry.dukhovskoy/topo/woa/'
fout = pthout+'etopo2woa025_bilinear.dat'
ijdm = idm*jdm
fid = open(fout,'rb')
try:
  fid.seek(0)
  HW = np.fromfile(fid,dtype='<f4', count=ijdm)
finally:
  fid.close()
HW = HW.reshape(jdm,idm)

# Find index to plot:
dmm = abs(lat-lat_plt)
jplt = np.argmin(dmm)
dmm = abs(lon-lon_plt)
iplt = np.argmin(dmm)

#dmm = S[:,jplt,iplt].squeeze()
#Swoa = dmm[np.where(~np.isnan(dmm))]
#dmm = T[:,jplt,iplt].squeeze()
#Twoa = dmm[np.where(~np.isnan(dmm))]
dmm = read_1prof(furlS,'s_an',jplt,iplt)
dmm[dmm>1.e10]=np.nan
Swoa = dmm[np.where(~np.isnan(dmm))]

dmm = read_1prof(furlT,'t_an',jplt,iplt)
dmm[dmm>1.e10]=np.nan
Twoa = dmm[np.where(~np.isnan(dmm))]

kobs = Twoa.shape[0]
Zwoa = ZZ[:kobs]


import mod_utils_wod as uwod
importlib.reload(uwod)
import mod_utils_rossby as ursb
importlib.reload(ursb)
import mod_calc_rrosby as clcr
importlib.reload(clcr)
import mod_solver as uslv

# ==========================================
#
#        START
#
# ==========================================

# Subset T/S profiles to specified depths to mimic
# observed values

Tobs = Twoa.copy()
Sobs = Swoa.copy()
Zobs = Zwoa.copy()
Tobs, Sobs = ursb.subset_TS2Z(Tobs,Sobs,Zobs,Zintrp)
Zobs = Zintrp

if np.min(Zobs) >= 0.:
  Zobs=-Zobs
#
# Find local depth
lon0 = lon[iplt]
lat0 = lat[jplt]
zbtm = HW[jplt,iplt] 

#
# QC check of T/S profiles
f_QC = ursb.run_QC(Zobs,Tobs,Sobs,zbtm)

# Vertical interp into Z, filling missing values
Tzz, Szz, ZZ = ursb.interp_TS_zlev(Zobs,Tobs,Sobs,zbtm,fintrp=1,add_srf=1)
#Tzz, Szz, ZZ = ursb.interp_TS_zlev(Zobs,Tobs,Sobs,zbtm,fintrp=0,add_srf=0)

f_plt = 0
if f_plt:
  import mod_plot_rsb as uplt
  importlib.reload(uplt)
  ctl = '{1},  Zbtm={0:.1f}'.format(zbtm,tfnm)
#      uplt.plot_prof(Zobs,Tobs,ctl=ctl) # original
  ctl2='T Interpolated'
  uplt.plot_2prof(Zobs,Tobs,ZZ,Tzz,ctl1=ctl,ctl2=ctl2)
  bottom_text(btx)

  ctl3 = '{1}, Zbtm={0:.1f}'.format(zbtm,tfnm)
  ctl4='S Interpolated'
  uplt.plot_2prof(Zobs,Sobs,ZZ,Szz,ctl1=ctl,ctl2=ctl4,fgn=2)
  bottom_text(btx)


# Compute N2 profile:
N2, Z_phi = clcr.calcN2_1D(Tzz,Szz,ZZ,lat0)

# Fill the gap from the deepest obs to the bottom
N20 = N2.copy()
Z_phi0 = Z_phi.copy()
ZZ0 = ZZ.copy()
N2, Z_phi, ZZ = clcr.intrpN2_bottom(N2,Z_phi,ZZ,zbtm)

# Compute Rossby R
Rrsb, Cphase, Phi = clcr.solve_SturmLiouv_1D(N2,Z_phi,ZZ,zbtm,lat0)

print('Rossby Radius = {0:.2f} km'.format(Rrsb))

f_testN = 1
if f_testN == 1:
  N2zF = uslv.runmn(N2,ZZ,mnwnd=7)  # Filter
  N_hrz = np.sqrt(N2)
  N_chr = 3600.*N_hrz   # cycles per hr
  NF_hrz = np.sqrt(N2zF)
  NF_chr = 3600.*NF_hrz

  plt.ion()
  fig4 = plt.figure(4,figsize=(8,8), constrained_layout=False)
  fig4.clf()
  ax1 = plt.axes([0.1, 0.1, 0.35, 0.8])
  plt.plot(N_chr,Z_phi)
  plt.plot(NF_chr,Z_phi)
  ax1.grid(True)
  ax1.set_xlabel('cycle/hr')
  ctl = 'N^2, Zbtm={0:.1f}, {1}'.format(zbtm,tfnm)
  ax1.set_title(ctl)

  ax2 = plt.axes([0.55, 0.1, 0.35, 0.8])
  ax2.plot(Phi,Z_phi)
  ax2.grid(True)
  ctl2 = '$\Phi$(z)'
  ax2.set_title(r'$\Phi$(z)')

  Txt = ['WOA as WOD cast']
  Txt.append('lat={0:.2f}N, lon={0:.2f}E'.format(lat0,lon0))
  Txt.append('Rossby(km)={0:.1f}'.format(Rrsb))
  Txt.append('PhaseSp(m/s)={0:.2f}'.format(Cphase))
  nt=len(Txt)

  x0 = 1
  y0 = 1
  ax3 = plt.axes([0.56,0.12,0.25,0.1])
  ax3.set(xlim=(x0,4.5),ylim=(y0,6))
  for itt in range(nt):
    ssl=Txt[itt]
    x1 = x0
    y1 = y0+itt
    ax3.text(x1,y1,ssl)
  ax3.axis('off')

  bottom_text(btx, pos=[0.01, 0.02])





