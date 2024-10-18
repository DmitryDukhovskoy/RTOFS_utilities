"""
  Use subset of SPEAR monthly climatology fields
  derived in derive_monthly_clim_spear.py

  Create daily fields from monthly by linear interpolation in time
  Note need to have +/- 1 month from the start/end of time for time interpolation 
"""

import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
#import pickle
import matplotlib.pyplot as plt
from yaml import safe_load

import mod_utils_ob as mutob
importlib.reload(mutob)


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
sys.path.append('./seasonal-workflow')
from boundary import Segment
import mod_time as mtime
import mod_mom6 as mmom6
import mod_utils as mutil

# Climatology derived for these years, started at mstart
# Inidicate start of the SPEAR forecast:
ens        = 1
yr_start   = 1993
mo_start   = 4
dnmb_start = mtime.datenum([yr_start,mo_start,1])
dv_start   = mtime.datevec(dnmb_start)

# Months of the f/cast for which daily data are being created
# Have -1 mont at the beginning and +1 mo at the end for interpolation
FMONTHS = [x for x in range(1,13)]
nmonths = len(FMONTHS)
TMM  = np.zeros((nmonths,4), dtype=int)  # f/cast time: year, month, Ndays in month
#dnmb = dnmb_start-15
dnmb = dnmb_start 
nmdays = 0
for imo in range(nmonths):
  dnmb = dnmb + nmdays
  dv = mtime.datevec(dnmb)
  nmdays = mtime.month_days(dv[1],dv[0])
  dnmb0 = int(mtime.datenum([dv[0],dv[1],15]))
  TMM[imo,:2] = dv[:2]
  TMM[imo,2] = int(nmdays)
  TMM[imo,3] = dnmb0
ndays = np.sum(TMM[:,2])

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

segments = [ Segment(1, 'north', hgrid, output_dir=outdir),
             Segment(2, 'east',  hgrid, output_dir=outdir),
             Segment(3, 'south', hgrid, output_dir=outdir),
             Segment(4, 'west',  hgrid, output_dir=outdir)]

nOB = len(segments)

# segments[0].__dict__

#static = xarray.open_dataset('/work/acr/spear/analysis/ocean_z.static.nc')
grid_spear = xarray.open_dataset('/work/Dmitry.Dukhovskoy/data/SPEAR/ocean_z.static.nc')
icegrid_spear = xarray.open_dataset('/work/Dmitry.Dukhovskoy/data/SPEAR/ice.static.nc')

fconfig = 'config_nep.yaml'
with open(fconfig) as ff:
  config = safe_load(ff)
spear_dir = os.path.join(config['filesystem']['spear_month_ens'], 'monthly_clim')

# Check if mapping indices exist, gmapi:
dirgmapi = config['filesystem']['spear_mom_gmapi']
flgmaph  = f'spear2mom_NEP_OB_gmapi_hpnt.nc'
flgmapu  = f'spear2mom_NEP_OB_gmapi_upnt.nc'
flgmapv  = f'spear2mom_NEP_OB_gmapi_vpnt.nc'
dflgmaph = os.path.join(dirgmapi, flgmaph)
dflgmapu = os.path.join(dirgmapi, flgmapu)
dflgmapv = os.path.join(dirgmapi, flgmapv)
# h-point indices
if not os.path.isfile(dflgmaph):
  print(f'Mapping indices hpnt are missing, {dflgmaph}, deriving ...')
  ds        = mutob.load_var(spear_dir, 'thetao', yr1=yr1, yr2=yr2, mstart=mstart)
  lon_spear, lat_spear = mutob.subset_spear_coord('h')
  mutob.spear2momOB_gmapi(segments, lon_spear, lat_spear, dflgmaph)
dsh = xarray.open_dataset(dflgmaph)

# u-point indices
if not os.path.isfile(dflgmapu):
  print(f'Mapping indices upnt are missing, {dflgmapu}, deriving ...')
  ds        = mutob.load_var(spear_dir, 'uo')
  lon_spear, lat_spear = mutob.subset_spear_coord('u')
  mutob.spear2momOB_gmapi(segments, lon_spear, lat_spear, dflgmapu)
dsu = xarray.open_dataset(dflgmapu)

# v-point indices
if not os.path.isfile(dflgmapv):
  print(f'Mapping indices vpnt are missing, {dflgmapv}, deriving ...')
  ds        = mutob.load_var(spear_dir, 'vo')
  lon_spear, lat_spear = mutob.subset_spear_coord('v')
  mutob.spear2momOB_gmapi(segments, lon_spear, lat_spear, dflgmapv)
dsv = xarray.open_dataset(dflgmapv)

icc = 0
dsetOB = xarray.Dataset()
for isgm in range(nOB):
  nsgm = isgm+1
  print(f'Processing lon/lat OB segment={nsgm}')
  dset   = mutob.derive_obsegm_lonlat(hgrid, segments, isgm)
  dsetOB = xarray.merge([dsetOB, dset])


# Year days for forecasts wrt to day 1 of the f/cast
# for monthly mean fields, days of the 15th day of the months
time_modays = np.zeros((12))
for imo in range(1,13):
  dnmb = TMM[imo-1,3]
  time_modays[imo-1] = dnmb - dnmb_start + 1 

# Time for daily fields wrt to day 1 of the f/cast:
time_days = np.array([x for x in range(1,ndays+1)])

# Check interpolated fields at a given location:
nsgm_chck = 2
ichck = 1201
kchck = 0
Tchck = np.zeros((2,12))  # for interpolation checking at 1 location
Tintp = np.zeros((2,ndays)) 
icc = -1
# Process by variable
# Prepare 1 year of data
npol = 3  # degree of interpolating polynom for temporal interp. monthly --> daily fields
spear_dir = config['filesystem']['nep_spear_subset'].\
                   format(year=dv_start[0], ens=ens)
for varnm in ['thetao', 'so']:
  # Load monthly fields for NEP subset SPEAR 
  flnm_spear = f'NEP_spear_{dv_start[0]}{dv_start[1]:02d}.{varnm}.nc'
  ds = mutob.read_spear_output(spear_dir, varnm, flnm_spear, fzint=True)

  # Spatial interpolation of 2D OB sections SPEAR --> NEP supergrid 
  for isgm in range(nOB):
    nsgm   = isgm+1
    print(f'\nProcessing {varnm} OB segment={nsgm}')
    INDX   = dsh[f'indx_segm{nsgm:03d}'].data
    JNDX   = dsh[f'jndx_segm{nsgm:03d}'].data
    # Interpolate to NEP MOM6 OB supergrid:
    dset   = mutob.derive_obsegm_3D(hgrid, ds, segments, isgm, varnm, 
                                    INDX, JNDX, time_steps=time_modays)
    # Time interpolation: monthly --> daily
    dsetI = mutob.interp_OBsegm_mnth2daily(dset, ndays, nsgm, varnm, npol=npol)
    dsetOB = xarray.merge([dsetOB, dsetI])

# ------------------------------------------
# Checking:
    icc += 1
    if nsgm == nsgm_chck:
      sfx   = f"segment_{nsgm:03d}"
      vards = f"{varnm}_{sfx}"

      dmm = dset[vards].data[:,kchck,:,:].squeeze()
      Tchck[icc,:] = dmm[:,ichck].squeeze()
      dmm = dsetI[vards].data[:,kchck,:,:].squeeze()
      Tintp[icc,:] = dmm[:,ichck].squeeze()
# ------------------------------------------

# ================================================
#  Plot time-interpolated fields for checking
# time series for selected locations
#  plot sections - end of the code
# ================================================
f_chckOB = False
if f_chckOB:
  time_intp = np.arange(1,ndays+1)
  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.55, 0.8, 0.35])
  ax1.plot(time_modays,Tchck[0,:],'.-')
  ax1.plot(time_intp,Tintp[0,:],'-')
  sttl = f'Daily interpolated thetao, P({npol})'
  ax1.set_title(sttl)
  ax1.grid('on')

  ax2 = plt.axes([0.1, 0.10, 0.8, 0.35])
  ax2.plot(time_modays,Tchck[1,:],'.-')
  ax2.plot(time_intp,Tintp[1,:],'-')
  sttl = f'Daily interpolated so, P({npol})'
  ax2.set_title(sttl)
  ax2.grid('on')


# Ssh - 1D sections
# Load ssh daily fields for NEP subset SPEAR 
varnm = 'ssh'
flnm_spear = f'NEP_spear_{dv_start[0]}{dv_start[1]:02d}.ssh_daily.nc'
ds = mutob.read_spear_output(spear_dir, varnm, flnm_spear, fzint=True)
for isgm in range(nOB):
  nsgm  = isgm+1
  print(f'Processing ssh OB segment={nsgm}')
  INDX   = dsh[f'indx_segm{nsgm:03d}'].data
  JNDX   = dsh[f'jndx_segm{nsgm:03d}'].data

  # Spatial interpolation from SPEAR --> NEP OB supergrid
  dset   = mutob.derive_obsegm_ssh(hgrid, ds, segments, isgm, INDX, JNDX, time_steps=time_days)
  dsetOB = xarray.merge([dsetOB, dset])


# Plot rot angles:
f_plt = False
if f_plt: mutob.plot_rotangle(grid_spear, icegrid_spear)
f_pltNEP = False
if f_pltNEP: mutob.plot_rotangleNEP(hgrid, hmask)

# UV fields
# Derive rotation angle, rad and components of the rotation matrix
# for rotating vectors onto true N/E grid from SPEAR
r2d = 180./np.pi
theta_rot, cosrot, sinrot = mutob.get_rotangle(icegrid_spear, fconfig, grid_spear)
print(f"Rotation angle for SPEAR min/max: {np.min(theta_rot)*r2d:6.2f} " + \
      f"/ {np.max(theta_rot)*r2d:6.2f}")

# Load monthly SPEAR data subset for NEP
flnmu_spear = f'NEP_spear_{dv_start[0]}{dv_start[1]:02d}.uo.nc'
flnmv_spear = f'NEP_spear_{dv_start[0]}{dv_start[1]:02d}.vo.nc'
ds_uo = mutob.read_spear_output(spear_dir, 'uo', flnmu_spear, fzint=True)
ds_vo = mutob.read_spear_output(spear_dir, 'vo', flnmv_spear, fzint=True)
for varnm in ['u', 'v']:
  for isgm in range(nOB):
    nsgm = isgm+1
    print(f'n\Processing {varnm} OB segment={nsgm}')

    if varnm == 'u':
      INDX  = dsu[f'indx_segm{nsgm:03d}'].data
      JNDX  = dsu[f'jndx_segm{nsgm:03d}'].data
    elif varnm == 'v':
      INDX  = dsv[f'indx_segm{nsgm:03d}'].data
      JNDX  = dsv[f'jndx_segm{nsgm:03d}'].data

    # Spatial interpolation SPEAR --> NEP OB supergrid
    dset   = mutob.derive_obsegm_uv(hgrid, ds_uo, ds_vo, segments, isgm, theta_rot,\
                                    varnm, INDX, JNDX, time_steps=time_modays)
    # Time interpolation: monthly --> daily
    dsetI = mutob.interp_OBsegm_mnth2daily(dset, ndays, nsgm, varnm, npol=npol)
    dsetOB = xarray.merge([dsetOB, dset])

for varnm in ['thetao','so','u','v']:
  for segm in [1,2,3,4]:
    vv  = f"{varnm}_segment_{segm:03d}"
    vdz = f"dz_{varnm}_segment_{segm:03d}"
    dm1 = f'lat_segment_{segm:03d}'
    dm2 = f'lon_segment_{segm:03d}'
    dsetOB[vv].attrs["coordinates"] = f"{dm1} {dm2}"
    dsetOB[vdz].attrs["coordinates"] = f"{dm1} {dm2}"

varnm='zos'
for segm in [1,2,3,4]:
  vv  = f"{varnm}_segment_{segm:03d}"
  dm1 = f'lat_segment_{segm:03d}'
  dm2 = f'lon_segment_{segm:03d}'
  dsetOB[vv].attrs["coordinates"] = f"{dm1} {dm2}"

dstart = f'{dv_start[0]}/{dv_start[1]:02d}/{dv_start[2]:02d}'
dsetOB.attrs["history"] = f"Created from SPEAR monthly T,S,U,V and daily SSH fields f/cast started {dstart} ens={ens:02d}"
dsetOB.attrs["code"] = f"/home/Dmitry.Dukhovskoy/python/setup_seasonal_NEP/create_daily_from_monthly_spear.py"

"""
  Add attributes for time var
"""
#dsetOB['time'] = np.arange(0, ntsteps, dtype='float')
dsetOB['time'].attrs['units'] = f'days since {dv_start[0]}-{dv_start[1]:02d}-{dv_start[2]:02d}'
dsetOB['time'].attrs['calendar'] = 'JULIAN'
dsetOB['time'].attrs['cartesian_axis'] = 'T'


encode = {
  'time': dict(dtype='float64', _FillValue=1.0e20)
}
for varnm in ['lon','lat']:
  for segm in [1, 2, 3, 4]:
    encode.update({
      f'{varnm}_segment_{segm:03d}': dict(dtype='float64', _FillValue=1.0e20)
    })
for varnm in ['zos']:
  for segm in [1, 2, 3, 4]:
    encode.update({
    f'{varnm}_segment_{segm:03d}': dict(_FillValue=1.0e20)
    })
for varnm in ['thetao', 'so', 'u', 'v']:
  for segm in [1, 2, 3, 4]:
    encode.update({
    f'{varnm}_segment_{segm:03d}': dict(_FillValue=1.0e20),
    f'dz_{varnm}_segment_{segm:03d}': dict(_FillValue=1.0e20)
    })

date_init = f'{dv_start[0]}{dv_start[1]:02d}{dv_start[2]:02d}'
pthoutp = gridfls['MOM6_NEP']['seasonal_fcst']['pthoutp']
fobc_out = os.path.join(pthoutp,f'OBCs_spear_daily_init{date_init}.nc')
print(f'Saving OBCs ---> {fobc_out}')
dsetOB.to_netcdf(fobc_out, 
                 format='NETCDF3_64BIT', 
                 engine='netcdf4', 
                 encoding=encode, 
                 unlimited_dims='time')

print('========  END ========')

# ===============================================================
f_chckOB = False
if f_chckOB: 
  varplt = 'thetao'
  nsgm = 2
  jday0 = 173   # day wrt to SPEAR forecast initial date
  dnmb0 = dnmb_start + jday0 - 1
  dv0   = mtime.datevec(dnmb0)
  # Find Start-end months for interpolated field
  kmin  = max(np.where(TMM[:,3] <= dnmb0)[0])
  kmax  = min(np.where(TMM[:,3] >= dnmb0)[0])
  if kmin == kmax: 
    if kmin > 0:
      kmin = kmin-1
    else:
      kmax = kmax+1
  
  dset_segm = segments[nsgm-1]
  xOB   = dset_segm.coords.lon.data
  yOB   = dset_segm.coords.lat.data
  if varplt == 'h':
    INDX   = dsh[f'indx_segm{nsgm:03d}'].data
    JNDX   = dsh[f'jndx_segm{nsgm:03d}'].data
  elif varplt == 'u':
    INDX   = dsu[f'indx_segm{nsgm:03d}'].data
    JNDX   = dsu[f'jndx_segm{nsgm:03d}'].data
  elif varplt == 'v':
    INDX   = dsv[f'indx_segm{nsgm:03d}'].data
    JNDX   = dsv[f'jndx_segm{nsgm:03d}'].data
  HHM = dstopo_nep['depth'].data
  HHM = np.where(HHM < 1.e-20, np.nan, HHM)
  HHM = -HHM
  HHM = np.where(np.isnan(HHM), 1., HHM)
  ds_topo_segm = mutob.segm_topo(nsgm, HHM, hgrid)
  Xsgm = ds_topo_segm['dist_supergrid'].data
  Xbtm = ds_topo_segm['dist_grid'].data
  Hbtm = ds_topo_segm['topo_segm'].data
  Hbtm = np.where(Hbtm > 0, 0., Hbtm)
  segm_nm = ds_topo_segm['segm_name'].data[0]
  Xbtm[-1] = Xsgm[-1]

  sfx   = f"segment_{nsgm:03d}"
  vards = f"{varplt}_{sfx}"
  dzvar = f"dz_{varplt}_{sfx}"

  flnm_spear = f'NEP_spear_{dv_start[0]}{dv_start[1]:02d}.{varplt}.nc'
  ds   = mutob.read_spear_output(spear_dir, varplt, flnm_spear, fzint=True)
  INDX   = dsh[f'indx_segm{nsgm:03d}'].data
  JNDX   = dsh[f'jndx_segm{nsgm:03d}'].data
  dset = mutob.derive_obsegm_3D(hgrid, ds, segments, nsgm-1, varplt,
                                    INDX, JNDX, time_steps=time_days)

  dZ     = dsetOB[dzvar].isel(time=jday0-1).data.squeeze()
  ZZ, ZM = mmom6.zz_zm_fromDZ(dZ)

  Fmo1 = dset[vards].isel(time=kmin).data.squeeze()
  Fmo2 = dset[vards].isel(time=kmax).data.squeeze()
  Fday = dsetOB[vards].isel(time=jday0-1).data.squeeze()

  clrmp, rmin, rmax, xl1, xl2 = mutob.OBsegm_clrmp_rminrmax(segm_nm, varplt, Xbtm=Xbtm) 
  btx = 'create_daily_from_monthly_spear.py'

  dstr = f'{TMM[kmin,0]}/{TMM[kmin,1]}'
  sttl = f'Monthly SPEAR {varplt} OB={segm_nm} {dstr}'
  mutob.plot_xsection(Fmo1, Xsgm, ZM, Hbtm, Xbtm, clrmp, rmin, rmax, xl1, xl2, sttl=sttl, fgnmb=1, btx=btx)

  dstr = f'{TMM[kmax,0]}/{TMM[kmax,1]}'
  sttl = f'Monthly SPEAR {varplt} OB={segm_nm} {dstr}'
  mutob.plot_xsection(Fmo2, Xsgm, ZM, Hbtm, Xbtm, clrmp, rmin, rmax, xl1, xl2, sttl=sttl, fgnmb=2, btx=btx)
    
  dstr = f'{dv0[0]}/{dv0[1]}/{dv0[2]}'
  sttl = f'Daily interp_P(2) {varplt} OB={segm_nm} {dstr}'
  mutob.plot_xsection(Fday, Xsgm, ZM, Hbtm, Xbtm, clrmp, rmin, rmax, xl1, xl2, sttl=sttl, fgnmb=3, btx=btx)


