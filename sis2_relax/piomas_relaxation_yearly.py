"""
  Create relax fields from PIOMAS monthly ice thickness and concentration

  usage: piomas_relaxation_yearly.py --yrs 1994 --yre 1995 --fsave 1

  monthly fields
  1901 - 2010
  https://psc.apl.uw.edu/research/projects/piomas-20c/

PIOMAS-20C is a sea ice thickness reconstruction covering the period 1901-2010. It is constructed using a coupled ice-ocean model using atmospheric forcing data from the ECMWF ERA-20C reanalysis to provide atmospheric forcing. Sea ice concentrations from the Hadley Center HadISST v2.0 data set are assimilated to constrain the model at the ice-edge. 

  All variables should have the following information for FMS subroutine 
  to process them correctly:

  axis information: src/mom6/src/MOM6/src/framework/MOM_io.F90:
  type :: axis_info
    character(len=32)  :: name = ""       !< The name of this axis for use in files
    character(len=256) :: longname = ""   !< A longer name describing this axis
    character(len=48)  :: units = ""      !< The units of the axis labels
    character(len=8)   :: cartesian = "N" !< A variable indicating which direction
                                          !! this axis corresponds with. Valid values
                                          !! include 'X', 'Y', 'Z', 'T', and 'N' for none.
    integer            :: sense = 0       !< This is 1 for axes whose values increase upward, or -1
                                          !! if they increase downward.  The default, 0, is ignored.
    integer            :: ax_size = 0     !< The number of elements in this axis
    real, allocatable, dimension(:) :: ax_data !< The values of the data on the axis [arbitrary]
  end type axis_info


"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
import matplotlib.pyplot as plt
import pickle
from yaml import safe_load
import argparse

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
import mod_time as mtime
import mod_mom6 as mmom6
import mod_utils as mutil
import mod_colormaps as mclrmps
import mod_misc1 as mmisc

parser = argparse.ArgumentParser()
parser.add_argument("--yrs",   help="year start to extract PIOMAS: 1993, ..., 2020", type=int)
parser.add_argument("--yre",   help="year end to extract PIOMAS: 1993, ..., 2020", type=int)
parser.add_argument("--fsave", help="flag > 0 to save the output", type=int)
args = parser.parse_args()

f_save = True
YRs = 1994
YRe = 1995    # make YRe=YRs to create 1 yr field with padded start/end of the year
file_type = 'monthly'  # monthly, daily, ... or clim
                       # for climatologies, do not need padded time - data will be recycled
                       # for monthly, daily, etc. need -dt and +dt at the beginn/end 


if args.yrs:
  YRs = args.yrs
if args.yre:
  YRe = args.yre
if args.fsave > 0:
  f_save = True
else:
  f_save = False


fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

# MOM6 NEP topo/grid:
run_name   = 'seasonal_fcst_daily'
pthtopo    = gridfls['MOM6_NEP'][run_name]['pthgrid']
fgrid      = gridfls['MOM6_NEP'][run_name]['fgrid']
ftopo_mom  = gridfls["MOM6_NEP"][run_name]["ftopo"]
outdir     = gridfls['MOM6_NEP'][run_name]['pthoutp']
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
# Hgrid lon. lat:
hlon, hlat  = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

pthsis  = gridfls['MOM6_NEP'][run_name]['pthsis']
pthdata = '/work/Dmitry.Dukhovskoy/data/PIOMAS_ice'
flthck = 'piomas20c.heff.1901.2010.v1.0.nc'
varthck = 'sit'
flconc  = 'piomas20c.area.1901.2010.v1.0.nc'
varconc = 'sic'

dflthkn = os.path.join(pthdata, flthck)
dflconc = os.path.join(pthdata, flconc)

ds_thkn = xarray.open_dataset(dflthkn)
LAT  = ds_thkn['Latitude'].data
LON  = ds_thkn['Longitude'].data

ds_conc = xarray.open_dataset(dflconc)

# Find gmapi - indices for bi-polar interpolation
import mod_regmom as mrmom
fgmapi  = f'PIOMAS_mom6_NEP_gmapi_{jdm}x{idm}.pkl'
dfgmapi = os.path.join(pthsis, fgmapi)

if os.path.isfile(dfgmapi):
  print(f'Loading gmapi <-- {dfgmapi}')
  with open(dfgmapi, 'rb') as fid:
   IMOM, JMOM, INDX, JNDX = pickle.load(fid) 
else:
  print('Searching gmapi for PIOMAS interpolation onto MOM6')
  jS = 570
  icc = -1
  IMOM = []
  JMOM = []
  for ii in range(idm):
    if ii%50 == 0:
      print(f' icc={icc} {ii/idm*100:.2f}% done ...')
    for jj in range(jS,jdm):
      if HH[jj,ii] >= 0:
        continue
      x0 = hlon[jj,ii]
      y0 = hlat[jj,ii]
      if y0 < 60.:
        continue
      if y0 < np.min(LAT) or y0 > np.max(LAT):
        continue

      icc += 1
      ixx, jxx = mrmom.find_gridpnts_box(x0, y0, LON, LAT, dhstep=1.)
      ixx = np.expand_dims(ixx, axis=0)
      jxx = np.expand_dims(jxx, axis=0)

      if icc == 0:
        INDX = ixx.copy()
        JNDX = jxx.copy()
      else:
        INDX = np.append(INDX, ixx, axis=0)
        JNDX = np.append(JNDX, jxx, axis=0)

      IMOM.append(ii)
      JMOM.append(jj)

  IMOM = np.array(IMOM)
  JMOM = np.array(JMOM)

  print(f'Saving gmapi --> {dfgmapi}')
  with open(dfgmapi, 'wb') as fid:
    pickle.dump([IMOM, JMOM, INDX, JNDX], fid)

# Create time array:
match file_type:
  case('clim'):
    TMPLT = np.zeros((12))
    for imo in range(1,13):
      dnmb0 = mtime.datenum([YRs,imo,15,12])
      TMPLT[imo-1] = dnmb0
  case('monthly'):
    TMPLT = []
    dnmb0 = mtime.datenum([YRs-1,12,15,12])
    TMPLT = [dnmb0]
    for YR in range(YRs,YRe+1):
      for imo in range(1,13):
        dnmb0 = mtime.datenum([YR,imo,15,12])
        TMPLT.append(dnmb0)
    dnmb0 = mtime.datenum([YRe+1,1,15,12])
    TMPLT.append(dnmb0)
    TMPLT = np.array(TMPLT)
  case _:
    raise Exception(f'relaxation input file for {file_type} has not been set up yet')

import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

nrec = len(TMPLT)
xdim = [x for x in range(0,idm)]
ydim = [x for x in range(0,jdm)]

H3d = np.zeros((nrec,jdm,idm))
C3d = np.zeros((nrec,jdm,idm))

#LMsk = HH.copy()
LMsk = np.where(HH<0, 1, 0)
icc  = -1
for dnmb in TMPLT:
  icc += 1

  dv = mtime.datevec(dnmb)
  yr0, mm0 = dv[:2]

  # Find record #: days since 1901-01-01 = day=1, index=0
  dnmb0 = mtime.datenum([yr0,mm0,1])
  dnmbR = mtime.datenum([1901,1,1])
  ndays = int(dnmb0-dnmbR) + 1
  #Time  = dset['time'].data  # np datetime array

  print(f'Processing {dv[0]}/{dv[1]}/{dv[2]}')
  Month = ds_thkn['month'].data
  Year  = ds_thkn['year'].data
  D     = np.sqrt((Month-mm0)**2 + (Year-yr0)**2)
  tindx = np.argmin(D)
  H2d   = ds_thkn[varthck].data[tindx,:].squeeze()  # thikness, m
  C2d   = ds_conc[varconc].data[tindx,:].squeeze()  #conc
  C2d   = np.where(C2d > 1., 1., C2d)

  # Get rid of nans - fill land:
  H2df = mmom6.fill_land3d(H2d)
  C2df = mmom6.fill_land3d(C2d)
  
  # interpolate onto MOM6 grid
  H2di = msisrlx.interp2Dfld(H2df, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)
  C2di = msisrlx.interp2Dfld(C2df, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)

  C2di = np.where(C2di>0.99, 1.0, C2di) # interp error 1--> 0.99, also caps max conc = 1

  H3d[icc,:,:] = H2di
  C3d[icc,:,:] = C2di

# Construct time array: days since the reference day:
dnmb_ref = mtime.datenum([1993,1,1])
dv_ref = mtime.datevec(dnmb_ref)
#dnmb0 = mtime.datenum([2003,12,15,12]) 
TM_ref = TMPLT - dnmb_ref 
if TM_ref[0] < 0:
  TM_ref[0] = 0.

# Construct data set:
dim1 = 'time'
dim2 = 'yh'
dim3 = 'xh'
ithknvar = 'ithkn'
iconcvar = 'iarea'  # partial area
darr_Hmom = xarray.DataArray(H3d, dims=(dim1, dim2, dim3), 
            coords={dim1: TM_ref, dim2: ydim, dim3: xdim})
darr_Cmom = xarray.DataArray(C3d, dims=(dim1, dim2, dim3), 
            coords={dim1: TM_ref, dim2: ydim, dim3: xdim})
dset_Hmom = xarray.Dataset({f'{ithknvar}': darr_Hmom})
dset_Cmom = xarray.Dataset({f'{iconcvar}': darr_Cmom})
dset_Ice = xarray.merge([dset_Hmom, dset_Cmom])

# Add attributes:
dset_Ice.attrs["history"] = f"Created from PIOMAS monthly ice fields {YRs}"
dset_Ice.attrs["code"] = "/home/Dmitry.Dukhovskoy/python/sis2_relax/piomas_relaxation_yearly.py"

dset_Ice[ithknvar].attrs["long_name"] = "Mean ice thickness"
dset_Ice[ithknvar].attrs["units"] = "meter"
dset_Ice[iconcvar].attrs["long_name"] = "Ice partial area, fraction"
dset_Ice[iconcvar].attrs["units"] = "unitless"

# Attributes for clim and time-varying fields:
match file_type:
  case('clim'):
    dset_Ice['time'].attrs['units'] = 'days since 0001-01-01'
    dset_Ice['time'].attrs['calendar'] = 'noleap'
    dset_Ice['time'].attrs['modulo'] = ' '
    dset_Ice['time'].attrs['cartesian_axis'] = 'T'
  case('monthly'):
    dset_Ice['time'].attrs['units'] = f'days since {dv_ref[0]}-{dv_ref[1]:02d}-{dv_ref[2]:02d}'
    dset_Ice['time'].attrs['calendar'] = 'gregorian'
    dset_Ice['time'].attrs['cartesian_axis'] = 'T'

dset_Ice['xh'].attrs['cartesian_axis'] = 'X'
dset_Ice['yh'].attrs['cartesian_axis'] = 'Y'

if f_save:
#  encoding = {rlx_name: {'_FillValue': None}}
  flout = f'PIOMAS_ithkn_iconc_{YRs}_{file_type}.nc'
  if not YRe == YRs:
    flout = f'PIOMAS_ithkn_iconc_{YRs}_{YRe}_{file_type}.nc'

  diclim = os.path.join(pthsis, flout)

  print(f'Saving PIOMAS climatology --> {diclim}')
  dset_Ice.to_netcdf(
       diclim,
       format='NETCDF3_64BIT',
       engine='netcdf4',
       unlimited_dims='time'
  )


check_rlx = False
if check_rlx:
  plt.ion()

  clrmp = mclrmps.colormap_temp2()
  rmin = 2.
  rmax = 20.

  # Stereographic Map projection:
  from mpl_toolkits.basemap import Basemap, cm
  m = Basemap(width=5000*1.e3,height=5000*1.e3, resolution='l',\
              projection='stere', lat_ts=50, lat_0=62, lon_0=-165)

  xR, yR = m(hlon, hlat)

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  m.drawcoastlines()
  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

  img = ax1.pcolormesh(xR, yR, RLXHR, cmap=clrmp, vmin=rmin, vmax=rmax)
#  img = ax1.pcolormesh(RLXHR, cmap=clrmp)

  ax1.set_title('Relaxation time, hrs')

  ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                     0.02, ax1.get_position().height])
  # extend: min, max, both
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
  ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.1f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  btx = 'relax_timescale.py'
  bottom_text(btx, pos=[0.2, 0.01])


