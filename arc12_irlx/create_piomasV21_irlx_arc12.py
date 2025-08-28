"""
  Create relaxation target fields from PIOMAS v21

  Need to have 4-vertices gmapi indices for bi-linear interpolation
  see: find_PIOMAS_to_ARC12_gmapi.py
  Find gmapi indices: 4 vertices of PIOMAS grid for each ARC12 grid
  for bilinear interpolation of PIOMAS --> ARC12
"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
import matplotlib.pyplot as plt
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
from mod_utils_fig import bottom_text
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

parser = argparse.ArgumentParser()
parser.add_argument("--yrs",   help="year start to extract PIOMAS: 1993, ..., 2020", type=int, required=True)
parser.add_argument("--yre",   help="year end to extract PIOMAS: 1993, ..., 2020", type=int)
parser.add_argument("--nyrs", help="number of years grouped in 1 relax. file, default=2", type=int)
parser.add_argument("--fsave", help="flag > 0 to save the output", type=int, required=True)
args = parser.parse_args()

file_type = 'monthly'  # monthly, daily, ... or clim
                       # for climatologies, do not need padded time - data will be recycled
                       # for monthly, daily, etc. need -dt and +dt at the beginn/end 

YRs  = args.yrs if args.yrs else None
nyrs = args.nyrs if args.nyrs else 2
YRe  = args.yre if args.yre else YRs
if args.fsave > 0:
  f_save = True
else:
  f_save = False


# ARC12 grid:
ptharc  = '/work/Dmitry.Dukhovskoy/ARC12/topo_grid'
dflarc  = os.path.join(ptharc,'ocean_hgrid.nc')
dfltopo = os.path.join(ptharc,'ocean_topog.nc')

ds_topo = xarray.open_dataset(dfltopo)
HH = -(ds_topo['depth'].data)
jdm, idm = HH.shape

assert HH[300,200] < 0., f'Check sign of topography, ocean pnts should be < 0'

hlon, hlat = mmom6.read_mom6grid(dflarc, grdpnt='hgrid') 

# PIOMAS LON/LAT:
import mod_regmom as mrmom
fgmapi  = f'PIOMAS_mom6_ARC12_gmapi_{jdm}x{idm}.npz'
pthgmapi = '/work/Dmitry.Dukhovskoy/ARC12/irlx'
pthsis = pthgmapi
dfgmapi = os.path.join(pthgmapi, fgmapi)

if os.path.isfile(dfgmapi):
  print(f'Loading gmapi <-- {dfgmapi}')
  data = np.load(dfgmapi)
  IMOM = data['IMOM']
  JMOM = data['JMOM']
  INDX = data['INDX']
  JNDX = data['JNDX']
else:
  raise Exception('Need to run find_PIOMAS_to_ARC12_gmapi.py to find gmapi first')

# Create time array:
def create_time_array(YR1, nyrs, file_type):
  YR2 = YR1+(nyrs-1)
  match file_type:
    case('clim'):
      TMPLT = np.zeros((12))
      for imo in range(1,13):
        dnmb0 = mtime.datenum([YR1,imo,15,12])
        TMPLT[imo-1] = dnmb0
    case('monthly'):
      TMPLT = []
      dnmb0 = mtime.datenum([YR1-1,12,15,12])
      TMPLT = [dnmb0]
      for YR in range(YR1,YR2+1):
        for imo in range(1,13):
          dnmb0 = mtime.datenum([YR,imo,15,12])
          TMPLT.append(dnmb0)
      dnmb0 = mtime.datenum([YR2+1,1,15,12])
      TMPLT.append(dnmb0)
      TMPLT = np.array(TMPLT)
    case _:
      raise Exception(f'relaxation input file for {file_type} has not been set up yet')

  return(TMPLT)

import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

pthdata = '/work/Dmitry.Dukhovskoy/data/PIOMAS_ice'
varthck = 'heff'
varconc = 'area'

for YR1 in range(YRs,YRe+1):
  YR2 = YR1 + (nyrs-1)
  print(f'Creating SIS2 relaxation for ARC12 {YR1}-{YR2}')
  TMPLT = create_time_array(YR1,nyrs,file_type)

  nrec = len(TMPLT)
  xdim = [x for x in range(0,idm)]
  ydim = [x for x in range(0,jdm)]

  H3d = np.zeros((nrec,jdm,idm))
  C3d = np.zeros((nrec,jdm,idm))

  #LMsk = HH.copy()
  LMsk = np.where(HH<0, 1, 0)
  icc  = -1
  YRold = 0
  for dnmb in TMPLT:
    icc += 1

    dv = mtime.datevec(dnmb)
    yr0, mm0 = dv[:2]

    if yr0 != YRold:
      flthck  = f'piomas_heff{yr0}_v21.nc'
      flconc  = f'piomas_area{yr0}_v21.nc'
      dflthkn = os.path.join(pthdata, flthck)
      dflconc = os.path.join(pthdata, flconc)

      ds_thkn = xarray.open_dataset(dflthkn)
      LAT  = ds_thkn['lat_scaler'].data
      LON  = ds_thkn['lon_scaler'].data

      ds_conc = xarray.open_dataset(dflconc)

      YRold = yr0

    # Find record #: monthly data
    #tindx = mm0-1

    Month = ds_thkn['month'].data
    Year  = ds_thkn['year'].data
    D     = np.sqrt((Month-mm0)**2 + (Year-yr0)**2)
    tindx = np.argmin(D)
    assert(D[tindx]==0), f"Requested {yr0}/{mm0} not found in {dflthkn}"
    print(f'Processing {dv[0]}/{dv[1]}/{dv[2]}, tindx={tindx}')
    H2d   = ds_thkn[varthck].data[tindx,:].squeeze()  # thikness, m
    C2d   = ds_conc[varconc].data[tindx,:].squeeze()  #conc
    #C2d   = np.where(C2d > 1., 1., C2d)

    # Get rid off the land values along the southern boundary:
    nbnd = 4
    for ik in range(nbnd):
      H2d[ik,:] = H2d[nbnd,:]
      C2d[ik,:] = C2d[nbnd,:]

    # Get rid of nans - fill land:
    H2df = mmom6.fill_land3d(H2d, land_mask=9999.9)
    C2df = mmom6.fill_land3d(C2d, land_mask=9999.9)
    C2df  = np.where(C2df > 1., 1., C2df)
    C2df  = np.where(C2df < 0.0, 0.0, C2df)

    # interpolate onto MOM6 grid
    H2di = msisrlx.interp2Dfld(H2df, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)
    C2di = msisrlx.interp2Dfld(C2df, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)

    # Fix the North Pole singularity that is not interpolated:
    iNP = 205
    jNP = 385
    if H2di[jNP,iNP] < 1.e-12:
      A = H2di[jNP-1:jNP+2,iNP-1:iNP+2]
      A = np.where(A<1.e-12, np.nan, A)
      H2di[jNP,iNP] = np.nanmean(A)
    if C2di[jNP,iNP] < 1.e-12:
      A = C2di[jNP-1:jNP+2,iNP-1:iNP+2]
      A = np.where(A<1.e-12, np.nan, A)
      C2di[jNP,iNP] = np.nanmean(A)

    # Fix the PIOMAS boundary swath of missing data near Greenland:
    iG1 = 270
    iG2 = 315
    jG1 = 608
    jG2 = 665
    for kfld in range(2):
      if kfld == 0:
        A2d = C2di.copy()
      else:
        A2d = H2di.copy()

      aa = A2d[jG1:jG2,iG1:iG2]
      aa = mmisc.box_fltr(aa)
      aa = mmisc.box_fltr(aa)
      A2d[jG1:jG2,iG1:iG2] = aa
      A2d = mmisc.box_fltr(A2d)
      A2d[HH>=0] = 0.

      if kfld == 0:
        C2di = A2d.copy()
      else:
        H2di = A2d.copy()

    C2di = np.where(C2di>0.99, 1.0, C2di) # interp error 1--> 0.99, also caps max conc = 1

    # Truncate very small ice thicknesses and concentrations
    iconc_min = 1.e-3
    ithkn_min = 1.e-2
    H2di[H2di<ithkn_min] = 0.
    C2di[C2di<iconc_min] = 0.

    # This should not happen but just in case:
    H2di[np.isnan(H2di)] = 0.
    C2di[np.isnan(C2di)] = 0.
     

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
  dset_Ice.attrs["history"] = f"Created from PIOMAS monthly ice fields {YR1}-{YR2}"
  dset_Ice.attrs["code"] = "/home/Dmitry.Dukhovskoy/python/sis2_relax/create_piomasV21_irlx_arc12.py"
  dset_Ice.attrs["info"] = f"Small ice conc < {iconc_min}, ice thkn < {ithkn_min} truncated to 0."

  dset_Ice[ithknvar].attrs["long_name"] = "Mean ice thickness or volume per m2"
  dset_Ice[ithknvar].attrs["units"] = "m3/m2"
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
    flout = f'PIOMASv21_ARC12_ithkn_iconc_{YR1}_{file_type}.nc'
    if not YR1 == YR2:
      flout = f'PIOMASv21_ARC12_ithkn_iconc_{YR1}_{YR2}_{file_type}.nc'

    diclim = os.path.join(pthsis, flout)

    print(f'Saving PIOMAS climatology --> {diclim}')
    dset_Ice.to_netcdf(
         diclim,
         format='NETCDF3_64BIT',
         engine='netcdf4',
         unlimited_dims='time'
    )


check_rlx = True
if check_rlx:
  plt.ion()

  varnm = 'ithkn'
  if varnm == 'iconc':
    clrmp = mclrmps.colormap_conc()
    clrmp.set_bad(color=[0.2, 0.2, 0.2])
    rmin = 0.
    rmax = 1.
    A2d = C2di
  elif varnm == 'ithkn':
    clrmp = mclrmps.colormap_ice_thkn()
    clrmp.set_bad(color=[0.2, 0.2, 0.2])
    rmin = 0.
    rmax = 5.
    A2d = H2di



  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])

  img = ax1.pcolormesh(A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.contour(HH,[0], linestyles='solid', colors=[(0,0,0)])
  ax1.contour(hlat,[40,50,60,64,70,80], linestyles='solid', colors=[(0.9,0.9,0.9)])

  ax1.set_title(f'PIOMAS interpolated to ARC12, {varnm}')
  ax1.axis('scaled')



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

  btx = 'create_piomasV21_irlx_arc12.py'
  bottom_text(btx, pos=[0.2, 0.01])

