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

# Start of the run - needed only for seasonal forecasts:
YRS    = 1993 # year start of the forecast
MMS    = 4
DDS    = 1
nens   = 2
dnmbR  = mtime.datenum([1994,3,3])  # day to plot

#dnmbS   = mtime.datenum([YRS,MOS,DDS])
#dayrun  = dnmbR - dnmbS + 1 # day to plot:
dvR     = mtime.datevec(dnmbR)

YR, MM, DD = dvR[:3]
jday  = int(mtime.date2jday(dvR[:3]))

expt    = "seasonal_fcst"
runname = f'NEPphys_frcst_climOB_{YRS}-{MMS:02d}-e{nens:02d}'
#expt    = 'NEP_BGCphys_GOFS'
#runname = 'NEP_physics_GOFS-IC'
#expt    = 'NEP_seasfcst_LZRESCALE'
#runname = 'NEPphys_LZRESCALE_climOB_1993_04-e02'


print(f'Expt: {expt} Run: {runname} Plot date: {dvR[0]}/{dvR[1]}/{dvR[2]}')
fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

if not expt == 'seasonal_fcst':
  YRS =  pthseas['MOM6_NEP'][expt]['year_start']
  MMS =  pthseas['MOM6_NEP'][expt]['month_start']
  DDS =  pthseas['MOM6_NEP'][expt]['day_start']

if expt == 'seasonal_fcst':
  pthfcst  = pthseas['MOM6_NEP'][expt]['pthoutp'].format(runname=runname)
else:
  dnmb0    = dnmbR
  dv0      = mtime.datevec(dnmb0)
  YR0, MM0, DD0 = dv0[:3]
  jday0    = int(mtime.date2jday([YR0,MM0,DD0]))
  pthfcst  = pthseas['MOM6_NEP'][expt]['pthoutp'].format(YY=YR0, MM=MM0)

pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
ndav       = pthseas['MOM6_NEP'][expt]['ndav']  # # of days output averaged
# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

# Find closest output:
if expt == 'NEP_BGCphys_GOFS':
  ocnfld = 'ocean'
else:
  ocnfld = 'oceanm'
if expt == 'seasonal_fcst':
  pthfcst = os.path.join(pthfcst,f'{ocnfld}_{dvR[0]}{dvR[1]:02d}')

if not os.path.isdir(pthfcst):
  print(f'not exist: {pthfcst}')

YR0, jday0, dnmb0, flname_out = manseas.find_closest_output(pthfcst, dnmbR, fld=ocnfld)
dv0  = mtime.datevec(dnmb0)
YR0, MM0, DD0 = dv0[:3]

# Averaging period:
dnmb_av1 = dnmb0 - np.floor(ndav/2)
#if dnmb_av1 < dnmbS: dnmb_av1=dnmbS
dnmb_av2 = dnmb_av1 + ndav-1

flocn_name = pthseas['MOM6_NEP'][expt]['focname'].format(YR=YR0, jday=jday0)
dfmom6 = os.path.join(pthfcst, flocn_name)

if not os.path.isfile(dfmom6):
  print(f'{dfmom6} not found ...')

dset = xarray.open_dataset(dfmom6)
SST  = dset['potT'].isel(time=0, zl=0).data

# Add observed ice conc:
# fields are downloaded from the Near-Real-Time NOAA/NSIDC 
# Climate Data Record of Passive Microwave Sea Ice Concentration 
# https://nsidc.org/data/g10016
f_cntrice = True
if f_cntrice:
  if expt == 'NEP_BGCphys_GOFS':
    outfld = 'ice'
  else:
    outfld = 'icem'

  if expt == 'seasonal_fcst':
    pthice = pthseas['MOM6_NEP'][expt]['pthoutp'].format(runname=runname)
    pthice = os.path.join(pthice,f'{outfld}_{dvR[0]}{dvR[1]:02d}')
  else:
    pthice  = pthseas['MOM6_NEP'][expt]['pthoutp'].format(YY=YR0, MM=MM0)

  flice_name = pthseas['MOM6_NEP'][expt]['ficename'].format(YR=YR0, jday=jday0)
  dfsis2 = os.path.join(pthice, flice_name)

  ds_ice   = xarray.open_dataset(dfsis2)
  I2d = ds_ice['siconc'].isel(time=0).data

  CNTR = manseas.derive_ice_contour(I2d, nmin=10)


dv_av1 = mtime.datevec(dnmb_av1)
yrs, mms, dds = dv_av1[:3]
dv_av2 = mtime.datevec(dnmb_av2)
yre, mme, dde = dv_av2[:3]

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

ss1 = f'{runname} SST \n'
sinfo = ss1 

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3300*1.e3, resolution='l',\
            projection='stere', lat_ts=55, lat_0=62, lon_0=-175)
xR, yR = m(hlon,hlat)


plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

img = m.pcolormesh(xR, yR, SST, cmap=clrmp, vmin=rmin, vmax=rmax)
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
      xc0, yc0 = mrtofs.interp_indx2lonlat(ii0, jj0, hlon, hlat)
      Xc[ipp] = xc0
      Yc[ipp] = yc0

    Xcm, Ycm = m(Xc,Yc)
    ax1.plot(Xcm, Ycm, linewidth=2, color=[1,0,0])


sttl = f"{runname} SST, ice edge avrg: {yrs}/{mms}/{dds}-{yre}/{mme}/{dde}"

ax1.set_title(sttl)

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

btx = 'plot_sst.py'
bottom_text(btx, pos=[0.2, 0.01])






