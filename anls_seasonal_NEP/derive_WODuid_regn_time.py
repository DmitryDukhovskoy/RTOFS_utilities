"""
  T/S Profile Observations - downloaded from WOD website

  https://www.ncei.noaa.gov/access/world-ocean-database-select/dbsearch.html
  Derive WOD UID for selected regions and time

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
import datetime
#import pickle
import xarray
import matplotlib.colors as colors
import matplotlib.mlab as mlab
#from netCDF4 import Dataset as ncFile
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
#import mod_mom6_valid as mom6vld
import mod_read_hycom as mhycom
import mod_misc1 as mmisc
import mod_time as mtime
import mod_mom6 as mmom6
import mod_WODdata as mwod
importlib.reload(mwod)


btx = 'derive_WODuid_regn_time.py'
regn = 'BeringSea'
# WOD data
pthwod = f'/work/Dmitry.Dukhovskoy/data/WOD_{regn}/'
pthoutp = '/work/Dmitry.Dukhovskoy/anls_output/NEPphys_frcst_dailyOB-expt02/data_misc/'

OBTYPES=['APB','CTD','PFL','OSD','XBT']

# Select lon, lat to search for WOD profiles
f_save = True
YR1    = 2015
YR2    = 2020
mo1    = 1
mo2    = 12

x0 = -170.
dx = 10.
y0 = 60.
dy = 10. 

print(f'Finding WOD UID for {regn} {YR1} {YR2} mo={mo1}-{mo2}\n')

Hmin = -200. # Depth limit, find sallow regions
icc = 0
for obtype in OBTYPES:
  prfx = 'ocldb1737657571.1064115.'
  furl = os.path.join(pthwod, obtype, f'{prfx}{obtype}.nc')
  Yobs, Xobs, TM, UID = mwod.search_UID(furl, x0, y0, dx, dy, \
                            YR1=YR1, YR2=YR2, mnth1=mo1, mnth2=mo2)

#DVctd = mtime.datevec1D(TMctd)
# Subset by depths, find observations in the regions shallower than Hmin
  for obtype in OBTYPES:
    Jctd, modify mwod.select_WODdepth to keep shallow regions 
  
#  Jctd, SIDctd = mwod.select_WODdepth(pthdata, UID, Hmin, LON, LAT, HH, qflag = True)

  nobs    = len(Xobs)
  darrx   = xarray.DataArray(Xobs, dims=(f"nobs_{obtype}"), \
                                   coords={f"nobs_{obtype}": np.arange(nobs)})
  darry   = xarray.DataArray(Yobs, dims=(f"nobs_{obtype}"), \
                                   coords={f"nobs_{obtype}": np.arange(nobs)})
  darrTM  = xarray.DataArray(TM,   dims=(f"nobs_{obtype}"), \
                                   coords={f"nobs_{obtype}": np.arange(nobs)})
  darrUID = xarray.DataArray(UID,  dims=(f"nobs_{obtype}"), \
                                   coords={f"nobs_{obtype}": np.arange(nobs)})

  if icc == 0:
    dset = xarray.Dataset({f"{obtype}_lon": darrx, f"{obtype}_lat": darry,\
                          f"{obtype}_TM": darrTM, f"{obtype}_UID": darrUID})
  else:
    dstmp = xarray.Dataset({f"{obtype}_lon": darrx, f"{obtype}_lat": darry,\
                           f"{obtype}_TM": darrTM, f"{obtype}_UID": darrUID})
    dset = xarray.merge([dset, dstmp])

  icc += 1

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

expt       = 'NEP_BGCphys_GOFS'
pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)

# Hgrid lon. lat:
LON, LAT = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

dflout = os.path.join(pthoutp, f'WOD_UID_{regn}_{YR1}_{YR2}.nc')

if f_save:
  print(f'Saving --> {dflout}')
  dset.to_netcdf(dflout,
      format='NETCDF3_64BIT',
      engine='netcdf4'
  )


f_pltobs = True
if not f_pltobs:
  return

# ===========================================
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.colors as colors
import matplotlib.mlab as mlab

CLRS = [[0.8, 0.4, 0],
        [0.2, 0.38, 1],
        [1, 0, 0.8],
        [0.5, 0.5, 0.5],
        [0, 0.6, 1],
        [0, 1., 0]]


# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3300*1.e3, resolution='l',\
            projection='stere', lat_ts=55, lat_0=62, lon_0=-175)

xR, yR = m(hlon, hlat)

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.05, 0.15, 0.8, 0.8])
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

hcntrs = [x for x in range(-200,0,20)]
m.contour(xR, yR, HH, hcntrs, colors=[(0.8,0.9,1)], linestyles='solid')

icc = 0
for obtype in OBTYPES:
  X = dset[f'{obtype}_lon'].data
  Y = dset[f'{obtype}_lat'].data  
  clr = CLRS[icc]
  Xm, Ym = m(X, Y)
  ax1.plot(Xm, Ym, 'o',  ms=3, markerfacecolor=clr, mec='none')
 
  icc += 1
 
sttl = f'WOD observationsi {regn} {YR1}-{YR2}'

ax2 = plt.axes([0.86,0.15,0.13,0.25])
icc = 0
y0 = 0.
dy = 0.1
for obtype in OBTYPES:
  clr = CLRS[icc]
  x0 = 0.1
  y0 = y0 + dy
  ax2.plot(x0, y0, 'o', ms=5, markerfacecolor=clr, mec='none')
  ax2.text(x0+dy, y0, obtype)
 
  icc += 1

ax2.set_xlim([x0-dy, x0+5*dy])
ax2.set_ylim([0.3*dy, y0+0.3*dy])
ax2.axis('off')




