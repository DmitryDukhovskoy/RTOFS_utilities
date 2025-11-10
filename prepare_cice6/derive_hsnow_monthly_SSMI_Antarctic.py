"""
  Derive monhtly clim of snow thickness in Antarctica 
  derived from AMSR 19 and 37 GHz microwave brightness temperatures
  https://earth.gsfc.nasa.gov/cryo/data/antarctic-snow-depth-sea-ice
  available data: 1993-2008
  daily data
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
import argparse

# Append custom module paths
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
#import mod_valid_utils as mvutil
import mod_colormaps as mclrmps

regn = 'south'  

parser = argparse.ArgumentParser()
parser.add_argument("--yrS", help="Start year to average, 1998,..., 2007", type=int, required=True)
parser.add_argument("--yrE", help="End year to average, default=2007", type=int)
args = parser.parse_args()

yrS    = args.yrS if args.yrS else None
yrE    = args.yrE if args.yrE else 2007

moS    = 1
moE    = 12
nmnths = moE-moS+1

assert moS <= moE, f'ERR: end month should be same or later than {moS}'

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


# Find year days for specified months:
def jdays(YR,moS,moE):
  jd1 = int(mtime.date2jday([YR,moS,1]))
  mo_days = mtime.month_days(moE,YR)
  jd2 = int(mtime.date2jday([YR,moE,mo_days]))

  return jd1, jd2

def read_hsnow(dflin):
  # Grid dimensions
  idim, jdim = 316, 332

  # Read the binary data
  with open(dflin, "rb") as fid:
    data = np.fromfile(fid, dtype=np.uint8)

  # Check size
  if data.size != jdim * idim:
    raise ValueError(f"Unexpected file size: expected {jdim * idim}, got {data.size}")

  A2D = data.reshape((jdim,idim))

  return A2D  

kyrs = 0
H3D = None
for YR in range(yrS,yrE+1):
  kyrs += 1

  match node_nm:
    case 'ppan':
      pthdata=f'/work/Dmitry.Dukhovskoy/data/snow_nasa/{YR}'
    case 'gaea' | 'dtn': 
      pthdata=f'/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/data/snow_nasa/{YR}'

  for MM in range (moS,moE+1):
    print(f'Processing {YR}/{MM}')

    jd1, jd2 = jdays(YR,MM,MM)
    icc = 0
    ASUM = None
    for jday in range(jd1,jd2+1):
      flnm = f's{YR}{jday:03d}.hs'
      dflin = os.path.join(pthdata,flnm)
      #print(f'Reading  {dflin}')

      # Some files may be missing, skip those
      if not os.path.isfile(dflin):
        print(f"File not found: {dflin}, skipping ...")
        continue  

      hs2D = read_hsnow(dflin).astype(float)
      # Flip the data to have correct orientation of Antarctica
      hs2D = np.flipud(hs2D)
      if ASUM is None:
        ASUM = hs2D.copy()
      else:
        ASUM = ASUM + hs2D
      icc += 1

    #ASUM = ASUM.astype(float)
    if icc > 1:
      HS2D = ASUM / float(icc)
    else:
      HS2D = ASUM.copy()

    # Land = 200
    HS2D[HS2D>190.]=np.nan   # thickns in cm !, 200 - land

    print(f'MM={MM}, Min/max snow (cm) = {np.nanmin(HS2D):.4f}/{np.nanmax(HS2D):.4f}')

    if H3D is None:
      jdim, idim = HS2D.shape
      H3D = np.zeros((nmnths,jdim,idim))

    # Sum of N years:
    H3D[MM-moS,:,:] = H3D[MM-1,:,:]+HS2D


if kyrs > 1:
  H3D = H3D/float(kyrs)

ntm, jdim, idim = H3D.shape


def read_NSIDC(YR,MM,DD,regn,pthnsidc,varnm):
  if regn == 'south':
    flnsidc = f"sic_pss25_{YR}{MM:02d}{DD:02d}_am2_v06r00.nc"
  else: 
    flnsidc = f"sic_psn25_{YR}{MM:02d}{DD:02d}_am2_v06r00.nc"
  
  with xarray.open_dataset(os.path.join(pthnsidc,flnsidc)) as ds_nsidc:
    A = ds_nsidc[varnm].data.squeeze()
    
  return A

# Get polar coordinates for Antarctica
# Same grid as in NSIDC ice conc. fields
# Southern Hemisphere Projection Based on WGS 1984
# https://nsidc.org/data/user-resources/help-center/guide-nsidcs-polar-stereographic-projection
YR = 2025
MM = 1
match node_nm:
  case 'ppan':
    pthnsidc = f'/work/Dmitry.Dukhovskoy/data/NRT_NOAA_NSIDC_seaconc/{YR}'
    pthout_snow = '/work/Dmitry.Dukhovskoy/data/snow_nasa/'
  case 'gaea' | 'dtn':
    pthnsidc = f'/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/data/NRT_NOAA_NSIDC_seaconc/{YR}'
    pthout_snow = '/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/data/snow_nasa/monthly_clim'

Xnsidc = read_NSIDC(YR,MM,1,regn,pthnsidc,'x')
Ynsidc0 = read_NSIDC(YR,MM,1,regn,pthnsidc,'y')
# NOAA NSIDC vertical axis is fliped upside down, so need to correct:
Ynsidc = np.flipud(Ynsidc0)

XX, YY = np.meshgrid(Xnsidc, Ynsidc, indexing='xy')
# Determine ellipsoid parameters from NSIDC information
# Note that Radius of ellipsoid WGS84 is typically referred to major semi-axis (equatorial radius)
if regn == 'south':
  slat0 = -70.  # latitude of 0 distortion, standard lat. 
  lon0  = -90.   # orientation of the 0 longitude wrt to X axis on polar grid
  flat_inv = 298.279411123064
  ax_maj = 6378273. 
  flat = 1./flat_inv  # flattening
  eccentr = np.sqrt(2*flat - flat**2)
  R_polar = ax_maj*(1.-flat)   # polar radius or semi-minor axis

if regn == 'south':
  LON, LAT = mmisc.convert_polarXY_lonlat(XX,YY, North=False, E=eccentr, RE=ax_maj, SLAT=slat0, LON0_dir=lon0)
  assert np.max(LAT) < 0., f"For southern hemisphere latitudes should be < 0"
  LON = -LON   # nor sure why but this makes sign of the longitudes right

# Double --> single precision
H3D = H3D.astype('float32')
LON = LON.astype('float32')
LAT = LAT.astype('float32')
Xnsidc = Xnsidc.astype('float32') 
Ynsidc = Ynsidc.astype('float32')

# Save to netCDF:
time_months = np.arange(1,13, dtype='int32')
darr_hs = xarray.DataArray(H3D, dims=("time","ypolar","xpolar"),\
                 coords={"time": time_months,\
                         "ypolar": Ynsidc,\
                         "xpolar": Xnsidc})
darr_xx = xarray.DataArray(Xnsidc, dims=("xpolar"),
                 coords={"xpolar": Xnsidc})
darr_yy = xarray.DataArray(Ynsidc, dims=("ypolar"),
                 coords={"ypolar": Ynsidc})
darr_lon = xarray.DataArray(LON, dims=("ypolar","xpolar"),
                 coords={"ypolar": Ynsidc,\
                         "xpolar": Xnsidc})
darr_lat = xarray.DataArray(LAT, dims=("ypolar","xpolar"),
                 coords={"ypolar": Ynsidc,\
                         "xpolar": Xnsidc})


dset_hs = xarray.Dataset(
    {
        "snow_depth": darr_hs,
        "lon": darr_lon,
        "lat": darr_lat,
        "xpolar": darr_xx,
        "ypolar": darr_yy,
    }
)

dset_hs['time'].attrs.update({
  "long_name": "months"
})
dset_hs['snow_depth'].attrs.update({
  "long_name": "snow depth on ice",
  "units": "cm",
  "land": "200"
})
dset_hs['xpolar'].attrs.update({
  "long_name": "polar stereographic X coordinates",
  "units": "m",
  "info": "Southern Hemisphere Projection Based on WGS 1984"
})
dset_hs['ypolar'].attrs.update({
  "long_name": "polar stereographic Y coordinates",
  "units": "m",
  "info": "Southern Hemisphere Projection Based on WGS 1984"
})
dset_hs['lon'].attrs.update({
  "long_name": "Longitudes",
  "units": "degrees",
})
dset_hs['lat'].attrs.update({
  "long_name": "Latitudes",
  "units": "degrees",
})

dset_hs.attrs.update({
  "title": f"NASA SSM/I-AMSR snow on ice monthly clim ({yrS}-{yrE})",
  "info": "Southern Hemisphere Snow depth files from SSM/I",
  "info2": "https://earth.gsfc.nasa.gov/cryo/data/antarctic-snow-depth-sea-ice",
  "institution": "NOAA NWS NCEP MDC",
  "source": "derive_hsnow_monthly_SSMI_Antarctic.py",
  "contact": "dmitry.dukhovskoy@noaa.gov",
  "region": regn,
  "Grid_idm_jdm": f"{idim}x{jdim}"
})

floutp = f"SSMI_Antarctic_hsnow_month_clim_{yrS}_{yrE}_{idim}x{jdim}.nc"
dflout = os.path.join(pthout_snow,floutp)
print(f"Saving --->   {dflout}")
dset_hs.to_netcdf(dflout, format="NETCDF4")


f_chck = False
if f_chck:
  plt.ion()

  fig1 = plt.figure(1,figsize=(12, 10))
  fig1.clf()  # Clear the figure

  units = 'cm'
  clrmp = mclrmps.colormap_temp()
  rmin = 0.
  #rmax = 20.
  rmax = 50.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])

  LON1 = LON.copy()
  LON2 = LON.copy()
  LON1 = np.where(LON1 < -175, np.nan, LON1)
  LON2 = np.where(LON2 > 172, np.nan, LON2)
  LON3 = np.where(LON < 0, LON+360., LON)
  LON3 = np.where(LON3 > 350., np.nan, LON3)
  lon_cntr1 = [x for x in range(-180,0,45)]  # grey -180:0
  lon_cntr2 = [x for x in range(45,178,45)]  # blue: 0 to 180 E
  lat_cntr = [x for x in range(-80,-20,10)]

  MM = 2
  sttl=f'hsnow MM={MM:02d}, {yrS}-{yrE}'

  hs2d = H3D[MM-1,:,:].squeeze()
  
  fig1 = plt.figure(1,figsize=(12, 10))
  fig1.clf()  # Clear the figure
  ax1 = plt.axes([0.08,0.1, 0.8,0.8])

  f_xy = True
  if f_xy: 
    img = ax1.pcolormesh(XX,YY,hs2d, cmap=clrmp, vmin=rmin, vmax=rmax)
    #ax1.invert_yaxis()
    # plot on XX,YY:
    cs = ax1.contour(XX,YY,LON1, lon_cntr1, linestyles='solid', colors=[(0.6,0.6,0.6)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs = ax1.contour(XX,YY,LON2, lon_cntr2, linestyles='solid', colors=[(0.,0.5,0.9)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs = ax1.contour(XX,YY,LAT, lat_cntr, linestyles='solid', colors=[(0.5,0.5,0.5)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs = ax1.contour(XX,YY,LON1,[0], linestyles='solid', colors=[(1,0.,0.)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs = ax1.contour(XX,YY,LAT,[-75], linestyles='solid', colors=[(1,0.,0.)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs = ax1.contour(XX,YY,LON3,[180], linestyles='solid', colors=[(0.6,0.,0.9)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    ax1.axis('scaled')

  else:
    img = ax1.pcolormesh(hs2d, cmap=clrmp, vmin=rmin, vmax=rmax)
    
    # Check longitudes:
    cs = ax1.contour(LON1, lon_cntr1, linestyles='solid', colors=[(0.6,0.6,0.6)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs2 = ax1.contour(LON2, lon_cntr2, linestyles='solid', colors=[(0.,0.5,0.9)], linewidths=1)
    ax1.clabel(cs2, inline=True, fontsize=10, fmt="%.1f")
    cs3 = ax1.contour(LON1,[0], linestyles='solid', colors=[(1,0.,0.)], linewidths=2)
    ax1.clabel(cs3, inline=True, fontsize=12, fmt="%.1f")
    ax1.contour(LON3,[180], linestyles='solid', colors=[(0.6,0.,0.9)], linewidths=1)
    cs = ax1.contour(LAT, lat_cntr, linestyles='solid', colors=[(0.5,0.5,0.5)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    cs = ax1.contour(LAT,[-75], linestyles='solid', colors=[(1,0.,0.)], linewidths=1)
    ax1.clabel(cs, inline=True, fontsize=10, fmt="%.1f")
    ax1.axis('scaled')

    #ax1.set_xlim([10,310])
    #ax1.set_ylim([20,330])
  ax1.set_title(sttl)

  # Colorbar
  # extend: min, max, both
  ax2 = fig1.add_axes([0.9,0.12,0.013,0.8])
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')

  ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  btx='derive_hsnow_monthly_SSMI_Antarctic.py'
  bottom_text(btx, pos=[0.1, 0.02])


