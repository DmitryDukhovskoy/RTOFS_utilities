"""
  OI SST high resolution fields
  https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html
# OpenDap to PSL data does not work on PPAN
# it works on Gaea
# I copied files PSL Linux
# [ddukhovskoy@linux256 noaa.oisst.v2.highres]$ pwd
#/Datasets/noaa.oisst.v2.highres
# ---> Niagara untrusted ---> Gaea --- gcp ---> PPAN
#
# Use subset data sets for NEP region
# see extract_OISST_NEPdomain.py
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

dnmb0 = mtime.datenum([1994,3,3])
dv0   = mtime.datevec(dnmb0)
YR, MM, DD = dv0[:3]
jday  = int(mtime.date2jday(dv0[:3]))


def read_field(furl,varnm):
  print("Reading {1} from {0}".format(furl,varnm))
  nc=ncFile(furl)
# lookup a variable
  dmm0 = nc.variables[varnm][:].data.squeeze()
  dmm = np.copy(dmm0)
  return dmm

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

# OpenDap to PSL data does not work on PPAN
#urlBase = 'http://psl.noaa.gov/thredds/dodsC/Datasets/noaa.oisst.v2.highres/'
#urlT    = f'sst.day.mean.{YR}.nc'
#furl = os.path.join(urlBase,urlT)
#T2d  = read_field(furl,'SST')
#ds = xarray.open_dataset(furl, chunks={})
#ds = xarray.open_dataset(furl)

pthsst = '/work/Dmitry.Dukhovskoy/data/OISST/'
#flsst  = os.path.join(pthsst,f'sst.day.mean.{YR}.nc')
flsst = os.path.join(pthsst,f'oisst_dayily_NEPsubset_{YR}.nc')
print(f'Reading {flsst}')
ds_sst = xarray.open_dataset(flsst)
lon1d  = ds_sst['lon'].data
lat1d  = ds_sst['lat'].data
SST    = ds_sst['sst'].isel(time=jday-1).data

LON, LAT = np.meshgrid(lon1d, lat1d)

# Add observed ice conc:
# fields are downloaded from the Near-Real-Time NOAA/NSIDC 
# Climate Data Record of Passive Microwave Sea Ice Concentration 
# https://nsidc.org/data/g10016
f_cntrice = True
if f_cntrice:
  pthnsidc = pthseas["NRT_NSIDC"]['pthdaily'].format(YR=YR)
  flnsidc  = f'seaice_conc_daily_nh_{YR}{MM:02d}{DD:02d}_f11_v04r00.nc'
  drflnsidc = os.path.join(pthnsidc, flnsidc)
  print(f'Reading NTR NSIDC ice conc: {drflnsidc}')
  dset_nsidc = xarray.open_dataset(drflnsidc)

  ICnrt = dset_nsidc['nsidc_nt_seaice_conc'].data[0,:].squeeze()
  Xnrt  = dset_nsidc['xgrid'].data
  Ynrt  = dset_nsidc['ygrid'].data

  ICnrt = np.where(ICnrt>1., np.nan, ICnrt)

# Flip NSIDC grid:
  ICnrt = np.flipud(ICnrt)
  Ynrt  = np.flipud(Ynrt)

# Convert Polar Coordinates to Geostatic coordinates (lon/lat)
  import mod_misc1 as mmisc
  XX, YY = np.meshgrid(Xnrt, Ynrt, indexing='xy')
  LONnrt, LATnrt = mmisc.convert_polarXY_lonlat(XX,YY)

# Make lon 0, 360 to match NEP grid
  LONnrt = np.where(LONnrt<0, LONnrt+360., LONnrt)

# Get ice edge contour in the Bering Sea
# get rid of ice in unneeded part of the domain
  ICnrt[:,150:] = np.nan
  ICnrt[:200,:] = np.nan

  CNTR = manseas.derive_ice_contour(ICnrt, nmin=5)

CLRS = [[0.8, 0.02, 0.6],
        [0.2, 0.38, 1],
        [1, 1, 1],
        [0., 0.8, 0.8],
        [0.4, 1, 0.7],
        [0., 0.8, 0],
        [0.3, 0.6, 0],
        [0.8, 0.8, 0],
        [1, 1, 0.5],
        [1, 0.9, 0.8],
        [1, 0.5, 0],
        [0.8, 0.4, 0.],
        [1, 0.65, 0.6],
        [1., 0., 0.],
        [0.7, 0.2, 0.1],
        [0.5, 0., 0.]]

import mod_colormaps as mclrmp
#rmin, rmax = mclrmp.minmax_clrmap(HpotZ, pmin=1, pmax=90)
rmin = -2.
rmax = 6.*abs(rmin)
#clrmp = mclrmp.colormap_temp2()
clrmp = mclrmp.colormap_posneg_uneven(CLRS)
clrmp.set_bad(color=[0.6,0.6,0.6])

ss1 = 'NOAA OI SST V2 High Resolution Dataset\n'
ss2 = f'NSIDC iconc: {drflnsidc}'
sinfo = ss1 + ss2

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3300*1.e3, resolution='l',\
            projection='stere', lat_ts=55, lat_0=62, lon_0=-175)
xR, yR = m(LON, LAT)


plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

img = m.pcolormesh(xR, yR, SST, cmap=clrmp, vmin=rmin, vmax=rmax)

sttl = f'OISST and ice edge from NSIDC ice conc {dv0[0]}/{dv0[1]}/{dv0[2]}'

ax1.set_title(sttl)

import mod_rtofs as mrtofs
if f_cntrice:
  ncc = len(CNTR)
  for icc in range(ncc):
    Ic = CNTR[icc][:,0]
    Jc = CNTR[icc][:,1]
# Get geodetic coordinates:
# use exact (float) indices to interpolate exact geodetic coordinate
    nic = len(Ic)
    Xc = np.zeros((nic))
    Yc = np.zeros((nic))
    for ipp in range(nic):
      ii0 = Ic[ipp]
      jj0 = Jc[ipp]
      xc0, yc0 = mrtofs.interp_indx2lonlat(ii0, jj0, LONnrt, LATnrt)
      Xc[ipp] = xc0
      Yc[ipp] = yc0

    Xcm, Ycm = m(Xc,Yc)
    ax1.plot(Xcm, Ycm, linewidth=2, color=[1,0,0])

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
# extend: min, max, both
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.1f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

ax3 = fig1.add_axes([0.02, 0.02, 0.8, 0.05])
ax3.text(0, 0, sinfo, fontsize=10)
ax3.axis('off')

btx = 'plot_OISST.py'
bottom_text(btx, pos=[0.2, 0.01])






