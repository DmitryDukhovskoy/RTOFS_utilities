"""
  Insert ice concentration (sea ice partial area)
  into CICE6 aice fields 
  Using NSIDC interpoalted fields

  See python/gfs_ice/interp_NSIDC_iconc_mesh025.py

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
from mpl_toolkits.basemap import Basemap, cm
import argparse

PPTHN = None
if 'PPTHN' not in locals() or PPTHN is None:
  cwd = os.getcwd()
  parts = cwd.split(os.sep)
  if 'python' in parts:
    idx = parts.index('python')
    PPTHN = os.sep + os.path.join(*parts[:idx + 1])
  else:
    raise RuntimeError("Directory 'python' not found in current working directory path.")

sys.path.extend([
    os.path.join(PPTHN, 'MyPython', 'hycom_utils'),
    os.path.join(PPTHN, 'MyPython', 'draw_map'),
    os.path.join(PPTHN, 'MyPython'),
    os.path.join(PPTHN, 'MyPython', 'mom6_utils')
])

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_cice6_utils as mc6util

rest_date = 20250103
rest_hr = 0
regn = 'south'

parser = argparse.ArgumentParser()
parser.add_argument("--rdate", help=f"restart date input file, default={rest_date}", type=int)
parser.add_argument("--rhr", help=f"input file, restart hour = 0, ..., 23, default={rest_hr}", type=int)
parser.add_argument("--rdate_out", help="output file, restart date if different from input", type=int)
parser.add_argument("--rhr_out", help="output file, restart hour if date is different from input", type=int)
parser.add_argument("--flrst_in", help="rest file in, otherwise name constructed from rest_date", type=str)
parser.add_argument("--flrst_out", help="new rest file, otherwise name constructed from rdate_out", type=str)
parser.add_argument("--regn", help=f"where icon incerted: south, north, global, default={regn}", type=str)
args = parser.parse_args()

rest_date = args.rdate if args.rdate else rest_date
rest_hr   = args.rhr if args.rhr else rest_hr
rest_date_out = args.rdate_out if args.rdate_out else rest_date
rest_hr_out   = args.rhr_out if args.rhr_out else rest_hr
flrst_in  = args.flrst_in if args.flrst_in else None
flrst_out = args.flrst_out if args.flrst_out else None

change_rest_time = (rest_date != rest_date_out) or (rest_hr != rest_hr_out)

# Get date numbers:
# Input restart file
dnmbR = mtime.rdate2datenum(rest_date*100+rest_hr)  # restart day nmb
yrR,mmR,ddR,hrR = mtime.datevec(dnmbR, round_hrs=True)[:4]
nsecR = hrR*3600

# Dates of the output fields in the new restart:
dnmbN = mtime.rdate2datenum(rest_date_out*100+rest_hr_out)
yrN, mmN, ddN, hrN = mtime.datevec(dnmbN, round_hrs=True)[:4]
nsecN = hrN*3600
 
syst_info = os.uname()
machine = syst_info.nodename
    
if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea compute node:", machine)
  node_nm = "gaea"
elif 'an' in machine:
  print("Running on PPAN node:", machine)  
  node_nm = "ppan"
else:
  print("Unknown machine:", machine)

fyaml = 'paths_ufs.yaml'
with open(fyaml) as ff:
  pths_ufs = safe_load(ff)

pthrest = os.path.join(pths_ufs[node_nm]["MOM6"]["pthrest"],'new')
pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]

# CICE parameters:
puny      = 1.e-11
c0        = 0.0
c1        = 1.0
c2        = 2.0
p5        = 0.5
Lsub      = 2.835e6    # latent heat sublimation fw (J/kg)
Lvap      = 2.501e6    # latent heat vaporization fw (J/kg)
Lfresh    = Lsub - Lvap # latent heat of melting of fresh ice (J/kg)
cp_ice    = 2106.       # specific heat of fresh ice (J/ kg/K)
rhos      = 330.        # density of snow (kg/m3)
hs_min    = 1.e-4       # min snow thickness for computing Tsno (m)
nsal      = 0.407
msal      = 0.573
min_salin = 0.1      # threshold for brine pocket treatment
saltmax   = 3.2        # max S at ice base
hg        = 1.e20    # bad values, land mask, etc.
nslyr     = 1   # snow layers
Tmin      = -100.   # minimum snow T


# Get MOM6 grid
pthgrid = pths_ufs[node_nm]["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

LON, LAT = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

# Interpolated NSIDC ice conc:
pthnsidc = os.path.join(pthdata,f"NRT_NOAA_NSIDC_seaconc/{yrN}")
if regn == 'south':
  RMsk = np.where(HH>=0, 0, 1)
  RMsk = np.where(LAT > -60., 0, RMsk)

print(f"old restart: {yrR}/{mmR:02d}/{ddR:02d}:{hrR:02d}")
print(f"new restart: {yrN}/{mmN:02d}/{ddN:02d}:{hrN:02d}")

fliceout = f'NSIDC_iconc_interp_mesh025_{jdm}x{idm}_{yrN}{mmN:02d}_{regn}.nc'
dfliceout = os.path.join(pthnsidc,fliceout)
print(f'Loading interpolated ice conc {dfliceout}')
with xarray.open_dataset(dfliceout) as dsint:
  AICEint = dsint['ice_conc'].isel(time=ddN-1).squeeze()

AICEint = np.where(RMsk == 0, np.nan, AICEint)


if flrst_in is None:
  flrst_in = f"cice_model.res.{yrR}{mmR:02d}{ddR:02d}.{nsecR:06d}.nc"
dflrst_in = os.path.join(pthrest, flrst_in)
print(f"Reading restart: {dflrst_in}")
ds_in = xarray.open_dataset(dflrst_in)
ds_out = ds_in.copy(deep=True)
ds_in.close()


# Insert aice and distribute/reduce proportionally over ice cats.:
assert nslyr == 1, f"Code needs to be modified for nslyr>1, nslyr={nslyr}"
aicen = ds_out['aicen'].data  # partial area by cats
vicen = ds_out['vicen'].data  # ice vol per m2 of grid cell by cats
vsnon = ds_out['vsnon'].data  # snow vol per m2 of ice area
qsnon = ds_out['qsno001'].data  # snow enthalpy by cats for 1 snow layer
Tsfcn = ds_out['Tsfcn'].data  # surface T in each cat.
ncat, jdim, idim = vsnon.shape

sice = {}
qice = {}
nilrs = 7   # ice layers
for i in range(1, nilrs+1):
  varnum = f"{i:03d}" 
  sice[varnum] = ds_out[f"sice{varnum}"].data
  qice[varnum] = ds_out[f"qice{varnum}"].data

# Aggregated ice partial area:
aice = np.sum(aicen, axis=0).squeeze()

# Select points to insert:
Jins, Iins = np.where((RMsk > 0) & (~np.isnan(AICEint)))
Xins = LON[Jins,Iins]
Yins = LAT[Jins,Iins]
npnts = len(Jins)

print(f"Found {npnts} points for insertion, min/max lat={np.min(Yins):.1f}/{np.max(Yins):.1f}"
       f" lon={np.min(Xins):.1f}/{np.max(Xins):.1f}")

def find_adj_icepnts(dlti, aice, i0, j0):
  jdm, idm = aice.shape
  i1 = max(i0 - dlti, 0)
  i2 = min(i0 + dlti, idm - 1)
  j1 = max(j0 - dlti, 0)
  j2 = min(j0 + dlti, jdm - 1)

  A = aice[j1:j2+1, i1:i2+1]
  jmm, imm = np.where(A > puny)

  # No grid points with ice found:
  if jmm.size == 0:
    return np.array([], dtype=int), np.array([], dtype=int)

  jice = jmm + j1
  iice = imm + i1

  return iice, jice

def find_adj_ocnpnts(dlti, aice, i0, j0):
  jdm, idm = aice.shape
  i1 = max(i0 - dlti, 0)
  i2 = min(i0 + dlti, idm - 1)
  j1 = max(j0 - dlti, 0)
  j2 = min(j0 + dlti, jdm - 1)

  A = aice[j1:j2+1, i1:i2+1]
  jmm, imm = np.where(A <= puny)

  # No grid points with no ice found:
  if jmm.size == 0:
    return np.array([], dtype=int), np.array([], dtype=int)

  jocn = jmm + j1
  iocn = imm + i1

  return iocn, jocn

# Note qsnon, qice < 0 !
vsnon_new = vsnon.astype(ds_out['vsnon'].dtype).copy()
qsnon_new = qsnon.astype(ds_out['qsno001'].dtype).copy()
aicen_new = aicen.astype(ds_out['aicen'].dtype).copy()
vicen_new = vicen.astype(ds_out['vicen'].dtype).copy()
Tsfcn_new = Tsfcn.astype(ds_out['Tsfcn'].dtype).copy()
qicen_new = {}
sicen_new = {}
for i in range(1, nilrs+1):
  varnum = f"{i:03d}"
  qicen_new[varnum] = qice[varnum].astype(ds_out[f"qice{varnum}"].dtype).copy()
  sicen_new[varnum] = sice[varnum].astype(ds_out[f"sice{varnum}"].dtype).copy()

Tfrz = -1.86243522  # ocean freez. T
Tsfc_max = -0.1  # max surf temp
dvol_sum = 0.
dlti = 2         # N of +/- i , j indices to search for ice pnts around i0,j0
hi_min = 0.05    # min ice thkn in each cat where aicen > 0 for noice --> ice case
hsnow_min = 0.01   # min snow thickness for noice --> ice case
print("Ice concentration insertion ...")
for ipp in range(npnts):
  if ipp%10000 == 0:
    print(f"   {ipp/npnts*100.:.2f}% done ...")
  j0 = Jins[ipp]
  i0 = Iins[ipp]

  # Note hsnow = vsn / aice for aice > 0
  # for cat n: vsn(n) = hsnow(n) * aice(n) 
  ai_old  = aice[j0,i0]       # aggreageted ice partial area 
  ain_old = aicen[:,j0,i0]    # partial areas by cats
  vsn_old = vsnon[:,j0,i0]    # snow volume per unit grid-cell area m2
  vin_old = vicen[:,j0,i0]    # ice volume per unit grid-cell area m2 
  #hin_old = np.divide(vin_old, ain_old, out=np.zeros_like(vin_old), where=ain_old != 0) # ice thkn or m3/m2_ice
  tsfcn_old = Tsfcn[:,j0,i0]  # surf T

  # Distribute new iconc proportionally by cats in snow vol m3/m2:
  ai_new = AICEint[j0,i0]
  if ai_new <= puny:
    ai_new = 0.
 
  if ai_new > 1.:
    ai_new = 1.

  # iconc change:
  if ai_new < puny:
    ain_new = ain_old * 0.
  else:
    if ai_old < puny:
      # all new iconc in cat 1:
      ain_new = ain_old * 0.
      ain_new[0] = ai_new
    else:
      dlt_ai = ai_new - ai_old
      wt = ain_old / np.sum(ain_old)
      dlt_ain = dlt_ai * wt
      ain_new = ain_old + dlt_ain

  assert np.sum(ain_new) <= 1., f"Check ain_new: sum>1: {np.sum(ain_new)}"

  aicen_new[:,j0,i0] = ain_new

  # Update ice enthalpy, Tsfcn, and salinity
  Tsf_mn = None
  iice = jice = iocn = jocn = None
  Vice_mn = None
  Vsn_mn  = None
  aice_case = None
  if ai_old <= puny and ai_new > puny:
    # Case: no ice --> ice, created ice in the grid cell
    aice_case = "noice2ice"

    # Update Tsfcn
    # Find N closest ice points:
    iice, jice = find_adj_icepnts(dlti, aice, i0, j0) 
    if len(iice) == 0:
      # no ice pnt adjacent to j0,i0:
      Tsf_mn = Tsfcn[:,j0,i0]*0.0 + Tsfc_max
      Vice_mn = hin_min * ain_old
    else: 
      # where ice - <= Tmax, where no ice = Tfrz
      Tsf_adj = Tsfcn[:,jice,iice]
      Tsf_mn  = np.nanmean(Tsf_adj, axis=1)
      Tsf_mn  = np.where(Tsf_mn > Tsfc_max, Tsfc_max, Tsf_mn)
      Tsf_mn  = np.where(ain_new < puny, Tfrz, Tsf_mn) 

      Vice_adj = vicen[:,jice,iice]
      Vice_mn = np.nanmean(Vice_adj, axis=1)
      hi_cell = np.max([np.sum(Vice_mn), hi_mn])    # m3/m2_cell or grid-cell mean ice thickness, m
      wt = ain_new / np.sum(ain_new)
      Vice_mn = hi_cell * wt

  elif ai_old > puny and ai_new > puny:
    # Case: ice --> updated ice conc
    aice_case = "ice2ice"
    Tsf_mn = np.where(tsfcn_old > Tsfc_max, Tsfc_max, tsfcn_old)
    Tsf_mn  = np.where(ain_new < puny, Tfrz, Tsf_mn)

    hi_cell = np.max([np.sum(vin_old), hi_min])    # m3/m2_cell or grid-cell mean ice thickness, m
    wt = vin_old / np.sum(vin_old)
    Vice_mn = hi_cell * wt

  elif ai_old > puny and ai_new < puny:
    # Case: ice --> no ice
    aice_case = "ice2noice"
    # Copy tsfc from adj grid cells
    iocn, jocn = find_adj_ocnpnts(dlti, aice, i0,j0)
    if len(iocn) == 0:
      # no ocn pnts:
      Tsf_mn = Tsfcn[:,j0,i0]*0.0 + Tfrz
    else:
      Tsf_adj = Tsfcn[:,jocn,iocn]
      Tsf_mn = np.nanmean(Tsf_adj, axis=1)
      Tsf_mn = np.where(Tsf_mn > Tsfc_max, Tsfc_max, Tsf_mn)
          
    Vice_mn = vin_old*0.0

  elif ai_old < puny and ai_new < puny:
    # Case: no ice --> no ice
    aice_case = "noice2noice"
    Tsf_mn = Tsfcn[:,j0,i0]
    Vice_mn = vin_old[:,j0,i0]*0.0
  
  else:
    raise Exception(f"Unexpected case for ai_old={ai_old} and ai_new={ai_new}")  

  # Should not happen but Checking if any NaN occur:
  if np.isnan(Tsf_mn).any():
    Tsf_mn = np.where(np.isnan(Tsf_mn), Tsfc_max, Tsf_mn)
  if np.isnan(Vice_mn).any():
    Vice_mn = np.where(np.isnan(Vice_mn), 0., Vice_mn)

  Tsfcn_new[:,j0,i0] = Tsf_mn
  vicen_new[:,j0,i0] = Vice_mn

  # Update sice00?, qice00?
  for ilr in range(1, nilrs+1):
    varnum = f"{ilr:03d}"
    # Ice salinity by layers - compute S profile using BZ99 formulation:
    sice_lr = mc6util.sice_lr_cice4(ilr, nilrs, ain_new) 
    sice_old = sice[varnum][:,j0,i0]
    sice_new = np.where(sice_old < puny, sice_lr, sice_old)
    sice_new = np.where(ain_new < puny, 0., sice_new)
    match aice_case:
      case "noice2ice":
        A = Stop

#------------
  # Update snow enthalpy: J/m3  
  # see icepack_therm_vertical.F90 in icepack
  #
  # snow enthalpy should be: qsn_min <= qsn <= qsn_max
  # In theory, qsn_max = -rhos_Lfresh (latent heat of metling at 0C)
  # Make it a little lower to keep snow from melting right away
  #hsn_new = vsn_new / ain
  qsn = qsnon[:,j0,i0]        # enthalpy, J/kg < 0
  qsn_min = -rhos * Lfresh + (Tmin + 0.01) * cp_ice * rhos  # enth. of the coldest possible snow
  qsn_max = -rhos * Lfresh - 0.01 * cp_ice * rhos  # a little colder than 0C snow
  qT0 = -Lfresh*rhos      # enth. of pure snow at 0C

  # Update new enthalpy of new snow:
  # Clip to min/max enthalpy, set to 0 where no snow:
  qsn_new = np.clip(qsn, qsn_min, qsn_max)
  qsn_new = np.where(ain <= puny, 0., qsn_new)      # no ice
  qsn_new = np.where(vsn_new <= 0, 0., qsn_new)     # no snow

  vtot_init = np.nansum(vsn)
  vtot_new  = np.nansum(vsn_new)
  #print(f"tot vsnon change = {vtot_new-vtot_init}") 

  dvol_sum = dvol_sum + (vtot_new-vtot_init)
  vsnon_new[:,j0,i0] = vsn_new
  qsnon_new[:,j0,i0] = qsn_new

  diff = np.nansum(vsn_new - vsnon[:, j0, i0])
  diff2 = np.nansum(vsn_new -vsn)
  diff3 = np.nansum(vsnon[:,j0,i0] - vsnon_new[:,j0,i0])
  if diff == 0 and abs(diff2) > 0:
    print(f"No change at {j0},{i0}, expected diff={diff2}")
  #else:
  #  print(f"diff={diff}, diff2={diff2}")  

  #assert abs(diff) > 0, f"no change at {j0},{i0}, expected diff={diff2}"

  if diff3 == 0 and abs(diff2) > 0:
    print(f"No change in the arrays at {j0},{i0}, expected diff={diff2}")

# Checking:
print(f"dvol_sum = {dvol_sum}")
total_vsnon_init = np.nansum(vsnon)
total_vsnon_new  = np.nansum(vsnon_new)
print(f"Tot snow volume change (m3/m2): {total_vsnon_new - total_vsnon_init}")

# Update data set:
#ds_out = ds_out.assign(vsnon=vsnon_new, qsno001=qsnon_new)
ds_out['vsnon'].values[:] = vsnon_new
ds_out['qsno001'].values[:] = qsnon_new


# Sanity checking:
assert "vsnon" in ds_out and "qsno001" in ds_out, "Missing updated snow fields vsnon and qsno001"
assert ds_out["vsnon"].shape == vsnon_new.shape, "Check shape of vsnon "
assert ds_out["qsno001"].shape == qsnon_new.shape, "Check shape of qsnon "

# Attributes:
from datetime import datetime
istep1_val = ds_out.attrs.get('istep1', None)
ds_out.attrs.update({
    "title": "CICE6 restart with inserted hsnow from SSM/I NASA gridded fields for S. Ocean",
    "source": "insert_hsnow_cice6_restart.py",
    "contact": "dmitry.dukhovskoy@noaa.gov",
    "istep1": np.int32(istep1_val) if istep1_val is not None else np.int32(0), 
    "myear": np.int32(yrN),
    "mmonth": np.int32(mmN),
    "mday": np.int32(ddN),
    "msec": np.int32(nsecN),
    "history": f"Modified {datetime.now().isoformat()}",
})

# Save:
if flrst_out is None:
  flrst_out = f"cice_model.res.{yrN}{mmN:02d}{ddN:02d}.{hrN:02d}.iconc.nc"
dflrst_out = os.path.join(pthrest,flrst_out)
print(f"Saving CICE restart --> {dflrst_out}")
ds_out.to_netcdf(dflrst_out, encoding={var: {'_FillValue': None} for var in ds_out.data_vars}, format='NETCDF3_64BIT')
ds_out.close()


f_plt = False
if f_plt:
  clrmp = mclrmps.colormap_conc()
  rmin = 0.
  rmax = 1.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  hlon = LON
  hlat = LAT   
 
  if regn == 'south':
    m = Basemap(projection='spstere',boundinglat=-50,lon_0=180,resolution='l')
  #lons, lats = m.makegrid(idim, jdim) # get lat/lons of ny by nx evenly spaced grid.
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)
  xl1 = -8.e6
  xl2 = -1.2e6
  yl1 = xl1
  yl2 = xl2

  xh, yh = m(hlon,hlat) # GFS coords

  plt.ion()
  fig1 = plt.figure(1, figsize=(8,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])

  m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
  m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
  img1 = ax1.pcolormesh(xh,yh,aice, cmap=clrmp, vmin=rmin, vmax=rmax)

  ax1.contour(xh,yh,HH,[0], linestyles='solid', colors=[(0.,0.,0.)], linewidths=1)

  ax1.set_xlim([xl1, xl2])
  ax1.set_ylim([yl1, yl2])
  ax1.invert_yaxis()
  ax1.invert_xaxis()

  # Plot pnt:
  x0 = hlon[j0,i0]
  y0 = hlat[j0,i0]
  xh0 = xh[j0,i0]
  yh0 = yh[j0,i0]
  ax1.plot(xh0,yh0,'o')

  # Colorbars
  ax3 = fig1.add_axes([0.2, 0.05, 0.6, 0.02])
  clb = plt.colorbar(img1, cax=ax3, orientation='horizontal', extend='max')
  ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax3.set_xticklabels(ax3.get_xticks())
  ticklabs = clb.ax.get_xticklabels()
  clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)










