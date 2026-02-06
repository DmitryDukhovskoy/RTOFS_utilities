"""
  Turning snwredist to 'ITDrdg' or 'ITDsd' requires additional 
  variables in the CICE6 restart file

  smice    , & ! tracer for mass of ice in snow (kg/m^3)
  smliq    , & ! tracer for mass of liquid in snow (kg/m^3)
  rsnw     , & ! snow grain radius (10^-6 m)
  rhos  snow density (kg/m^3)

  Add missing variables to CICE6 restart

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
importlib.reload(mc6util)


def main():
  flrst_in = 'cice_model.res.20250103.00.iconc_thkn.snow.nc'
  sfx_end = 'snphys'
  yrR = mmR = ddR = hrR = None
  yrN = mmN = ddN = hrN = None

  parser = argparse.ArgumentParser()
  parser.add_argument("--pth_in", help="input restart directory with original file", type=str, required=True)
  parser.add_argument("--flrst_in", help=f"rest file in", required=True, type=str)
  parser.add_argument("--pth_out", help="output restart directory where new file be dumped", 
                      type=str, required=True)
  parser.add_argument("--flrst_out", help="new rest file, otherwise flrst_in + snowphys", type=str)
  args = parser.parse_args()

  flrst_in   = args.flrst_in  if args.flrst_in  is not None else flrst_in
  flrst_out  = args.flrst_out if args.flrst_out is not None else None
  pthrst_in  = args.pth_in    if args.pth_in    is not None else None
  pthrst_out = args.pth_out   if args.pth_out   is not None else None

  # Derive restart dates assuming file nameing is cice_restart.res.YYYYMMDD.XX[XXX]
  # or YYYYMMDD.<time>.---.nc
  yrR, mmR, ddR, hrR, mintR = mc6util.get_date_filename(flrst_in)
  rest_date = int(yrR*1e4 + mmR*100 + ddR)
  rest_hr = hrR

  # Output restart - same date
  yrN, mmN, ddN, hrN, mintN = yrR, mmR, ddR, hrR, mintR
  rest_date_out = int(yrN*1e4 + mmN*100 + ddN)
  rest_hr_out = hrN

  print(f"Restart date input:  {rest_date}:{rest_hr}")
  print(f"Restart date output: {rest_date_out}:{rest_hr_out}")

  dnmbR = mtime.datenum([yrR,mmR,ddR,hrR])
  dnmbN = dnmbR

  nsecR = hrR*3600
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

  # Default parameters (icepack_parameters.F90)
  nslyr      = 1          # snow layers
  rhos       = 330.       # density of snow (kg/m3)
  puny       = 1.e-11
  c0         = 0.0
  c1         = 1.0
  c2         = 2.0
  p5         = 0.5
  Lsub       = 2.835e6    # latent heat sublimation fw (J/kg)
  Lvap       = 2.501e6    # latent heat vaporization fw (J/kg)
  Lfresh     = Lsub - Lvap # latent heat of melting of fresh ice (J/kg)
  cp_ice     = 2106.       # specific heat of fresh ice (J/ kg/K)
  rsnw_fall  = 54.526     # radius of new snow (10^-6 m)
  rsnw_tmax  = 1500.0     # maximum snow radius (10^-6 m)
  rhosnew    =  100.0     # new snow density (kg/m^3)
  rhosmin    =  100.0     # minimum snow density (kg/m^3)
  rhosmax    =  450.0     # maximum snow density (kg/m^3)
  windmin    =   10.0     # minimum wind speed to compact snow (m/s)
  drhosdwind =   27.3     # wind compaction factor for snow (kg s/m^4)
  snwlvlfac  =    0.3     # fractional increase in snow depth for bulk redistribution
  snw_growth_wet = 4.22e5 # wet metamorphism parameter (um^3/s) 1.e18 * 4.22e-13 (Oleson 2010)
  drsnw_min  =    0.0     # minimum snow grain growth factor
  snwliq_max =    0.033   # irreducible saturation fraction
                          #   0.033 (Anderson 1976)
                          # 0.09 to 0.1  (Denoth et al, 1979 & Brun 1989)

  print(f"old restart: {yrR}/{mmR:02d}/{ddR:02d}:{hrR:02d}")
  print(f"new restart: {yrN}/{mmN:02d}/{ddN:02d}:{hrN:02d}")

  #if flrst_in is None:
  #  flrst_in = f"cice_model.res.{yrR}{mmR:02d}{ddR:02d}.{nsecR:06d}.nc"

  dflrst_in = os.path.join(pthrst_in, flrst_in)
  print(f"Reading restart: {dflrst_in}")
  ds_in = xarray.open_dataset(dflrst_in)
  ds_out = ds_in.copy(deep=True)
  ds_in.close()

  # Add new variables:
  dims = ds_out['qsno001'].dims
  coords = ds_out['qsno001'].coords
  #attrs  = ds_out['qsno001'].attrs  

  assert nslyr == 1, f"Code needs to be modified for nslyr>1, nslyr={nslyr}"
  aicen = ds_out['aicen'].data  # partial area by cats
  vicen = ds_out['vicen'].data  # ice vol per m2 of grid cell by cats
  vsnon = ds_out['vsnon'].data  # snow vol per m2 of ice area
  qsnon = ds_out['qsno001'].data  # snow enthalpy by cats for 1 snow layer
  Tsfcn = ds_out['Tsfcn'].data  # surface T in each cat.
  ncat, jdim, idim = vsnon.shape

  # New fields:
  dtype = ds_out['qsno001'].dtype
  smliq = np.zeros_like(qsnon, dtype=dtype)
  smice = np.zeros_like(qsnon, dtype=dtype)
  rhosn = np.zeros_like(qsnon, dtype=dtype)
  rsnw  = np.zeros_like(qsnon, dtype=dtype)

  # Typical snow-to-liquid ratios by snow temperature:
  # snow:liquid = 10:1 for ~0 snow, and ~0 for T<0
  # For CICE6 snow, asumme no liquid fraction if T snow < Tliq
  # snow grain radius 
  # Typical effective snow grain radii in the Arctic 
  # vary widely, from around 50-100 µm (micrometers) 
  # for fresh snow to over 1000 µm (1 mm) for old, wet, or melting snow
  # Dang et al., JGR Atmos, 2017, "Measurements of light-absorbing particles in snow ..."

  Tliq = -0.05
  fliq_cold = 0.  # Fraction of liquid part in snow for cold T
  fliq_warm = 0.02 # Fraction of liquid part in snow, keep low to prevent rapid refreezing
  rsnw_warm = 500.     # grain size micro-m for warm snow
  rsnw_cold = 100.

  # Constant snow density is assumed 
  # which is true for restart for the no snow physics/metamorphosis cases
  for kk in range(ncat):
    aik  = aicen[kk,:,:]
    qsk  = qsnon[kk,:,:]
    tsk  = (qsk + rhos*Lfresh) / (cp_ice*rhos)  # should be close to surface T for 1lr snow
    warm = tsk > Tliq
    valid = aik > puny

    smliq[kk, :, :] = np.where(warm & valid, rhos * fliq_warm, 0.0)
    smice[kk, :, :] = np.where(valid, rhos - smliq[kk, :, :], 0.0) # to guarantee smliq+smice = rhos
    rhosn[kk, :, :] = np.where(valid, rhos, 0.0)
    rsnw[kk, :, :] = np.where(
        valid,
        np.where(warm, rsnw_warm, rsnw_cold),
        0.0
    )


  # Add new variables:
  ds_out['smliq001'] = xarray.DataArray(smliq, dims=dims, coords=coords)
  ds_out['smice001'] = xarray.DataArray(smice, dims=dims, coords=coords)
  ds_out['rhos001']  = xarray.DataArray(rhosn, dims=dims, coords=coords)
  ds_out['rsnw001']  = xarray.DataArray(rsnw, dims=dims, coords=coords)

  for v in ['smliq001', 'smice001', 'rhos001', 'rsnw001']:
      ds_out[v].encoding = ds_out['qsno001'].encoding.copy()

  # Physical consistency
  np.testing.assert_allclose(
      ds_out['smliq001'] + ds_out['smice001'],
      ds_out['rhos001']
  )

  print("Created snow fields:")
  for k in range(1,ncat+1):
    a2d = smliq[k-1,:,:]
    print(f"cat={k} smliq001 min/max: {np.nanmin(a2d)}/{np.nanmax(a2d)}")
    a2d = smice[k-1,:,:]
    print(f"cat={k} smice001 min/max: {np.nanmin(a2d)}/{np.nanmax(a2d)}")  
    a2d = rhosn[k-1,:,:]
    print(f"cat={k} rhos001  min/max: {np.nanmin(a2d)}/{np.nanmax(a2d)}")
    a2d = rsnw[k-1,:,:]
    print(f"cat={k} rsnw001  min/max: {np.nanmin(a2d)}/{np.nanmax(a2d)}\n")

  # Attributes:
  from datetime import datetime
  istep1_val = ds_out.attrs.get('istep1', None)
  ds_out.attrs.update({
      "title": f"CICE6 restart added snow parameters fields for snow physics",
      "info": f"Modified restart: {flrst_in}",
      "source": "add_snowvar_snowphys_cice6rest_global.py",
      "istep1": np.int32(istep1_val) if istep1_val is not None else np.int32(0), 
      "myear": np.int32(yrN),
      "mmonth": np.int32(mmN),
      "mday": np.int32(ddN),
      "msec": np.int32(nsecN),
      "history3": f"Modified {flrst_in} {datetime.now().isoformat()}",
  })

  # Save:
  if flrst_out is None:
    base, ext = os.path.splitext(flrst_in)
    flrst_out = f"{base}.{sfx_end}{ext}"

  dflrst_out = os.path.join(pthrst_out, flrst_out)
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
    xl1 = -9e6
    xl2 = -0.8e6
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

if __name__ == "__main__":
  main()

