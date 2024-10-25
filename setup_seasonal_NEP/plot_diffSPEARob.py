"""
  Plot difference of NEP OBs derived from  subset of SPEAR monthly fields
  derived in derive_monthly_clim_spear.py
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
import mod_utils as mutil
from boundary import Segment
import mod_time as mtime
import mod_mom6 as mmom6
from mod_utils_fig import bottom_text

# Climatology derived for these years, started at mstart
extract_new = False  # extract new set of SPEAR data, otherwise load already extracted 
yr_start = 1993
mo_start = 4
ens      = 1
mo_plt   = 1  # f/cast month to plot wrt to start time, i.e. from 1 to 12
varplt   = 'thetao'  # var to plot
dnmb_start = mtime.datenum([yr_start,mo_start,1])
dv_start   = mtime.datevec(dnmb_start)
run_name   = 'seasonal_fcst_daily'

ENSR = [1,2,3,4,5,6] # ensemble runs
nens = len(ENSR)

#indir = Path('/work/acr/spear/processed/ensmean')
#outdir = Path('/work/acr/spear/climatology')

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
pthtopo    = gridfls['MOM6_NEP'][run_name]['pthgrid']
fgrid      = gridfls['MOM6_NEP'][run_name]['fgrid']
ftopo_mom  = gridfls["MOM6_NEP"][run_name]["ftopo"]
outdir     = gridfls['MOM6_NEP'][run_name]['pthoutp']
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc')) 
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom  = os.path.join(pthtopo, fgrid)
# Hgrid lon. lat:
hlon, hlat  = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

segments = [ Segment(1, 'north', hgrid, output_dir=outdir),
             Segment(2, 'east',  hgrid, output_dir=outdir),
             Segment(3, 'south', hgrid, output_dir=outdir),
             Segment(4, 'west',  hgrid, output_dir=outdir)]

nOB = len(segments)

date_init = f'{dv_start[0]}{dv_start[1]:02d}{dv_start[2]:02d}'
pthoutp = gridfls['MOM6_NEP'][run_name]['pthoutp']
fobc_out = os.path.join(pthoutp,f'OBCs_spear_mnth_check{date_init}.nc')
if extract_new:
  dsetOB = mutob.monthlyOB_from_SPEAR_ensmbls(dnmb_start, varplt, run_name, ENSR )
  print(f'Saving SPEAR OBCs from {nens} ensemlbe runs ---> {fobc_out}')
  dsetOB.to_netcdf(fobc_out,
                 format='NETCDF3_64BIT',
                 engine='netcdf4',
                 unlimited_dims='time')
else:
  print(f'Loading SPEAR OBCs from {nens} ensemlbe runs ---> {fobc_out}')
  dsetOB = xarray.open_dataset(fobc_out)

# ===================
# Plotting
# ===================

dltx = 0.06
dlty = 0.07
dx = 0.27
dy = 0.27
xl = 0.06
yb = 0.15
FPOS = [[xl, yb+dy+dlty, dx, dy],
        [xl+dx+dltx, yb+dy+dlty, dx, dy],
        [xl+2*(dx+dltx), yb+dy+dlty, dx, dy],
        [xl, yb, dx, dy],
        [xl+dx+dltx, yb, dx, dy],
        [xl+2*(dx+dltx), yb, dx, dy]]



from matplotlib.patches import Polygon

plt.ion()

dnmb0 = TMM[mo_plt-1,3]
dv0   = mtime.datevec(dnmb0)
HHM = dstopo_nep['depth'].data
HHM = np.where(HHM < 1.e-20, np.nan, HHM)
HHM = -HHM
HHM = np.where(np.isnan(HHM), 1., HHM)


iref = 0   # Reference ensemble for computing differences
Nsegms = len(segments)
for isgm in range(Nsegms):
  nsgm  = isgm+1
  fgnmb = isgm+1
  print(f'Plotting segment {nsgm}')
  fig1 = plt.figure(fgnmb,figsize=(12,9))
  plt.clf()

  ds_topo_segm = mutob.segm_topo(nsgm, HHM, hgrid)
  Xsgm = ds_topo_segm['dist_supergrid'].data
  Xbtm = ds_topo_segm['dist_grid'].data
  Hbtm = ds_topo_segm['topo_segm'].data
  Hbtm = np.where(Hbtm > 0, 0., Hbtm)
  segm_nm = ds_topo_segm['segm_name'].data[0]
  Xbtm[-1] = Xsgm[-1]

  # Reference field:
  sfx   = f"e{iref+1:02d}_segment_{nsgm:03d}"
  vards = f"{varplt}_{sfx}"
  dzvar = f"dz_{varplt}_segment_{nsgm:03d}"
  F2dR  = dsetOB[vards].isel(time=mo_plt-1).data.squeeze()
  
  
  for ie in range(nens):
    ax1 = plt.axes(FPOS[ie])
    sfx   = f"e{ie+1:02d}_segment_{nsgm:03d}"
    vards = f"{varplt}_{sfx}"
    dzvar = f"dz_{varplt}_segment_{nsgm:03d}"
    dZ     = dsetOB[dzvar].isel(time=mo_plt-1).data.squeeze()
    ZZ, ZM = mmom6.zz_zm_fromDZ(dZ)
    F2d = dsetOB[vards].isel(time=mo_plt-1).data.squeeze()
    dF  = F2d - F2dR
   
    clrmp, rmin, rmax, xl1, xl2 = mutob.OBsegm_clrmp_rminrmax(segm_nm, varplt, Xbtm=Xbtm)
    clrmp  = mutil.colormap_ssh(nclrs=200)
    rmin = -1.
    rmax = 1.
    cntr1 = [x/10 for x in range(-10, 10, 1)]

    btx = 'plot_diffSPEARob.py'

    dstr = f'{TMM[mo_plt-1,0]}/{TMM[mo_plt-1,1]}'
    dstart = f'{dv_start[0]}/{dv_start[1]}/{dv_start[2]}'
    sttl = f'ens={ie+1:02d}'

    if ie == 0:
      btx0     = btx
      stxt     = f'Monthly SPEAR init: {dstart} diff {varplt} OB={segm_nm} wrt OB={iref+1}, Plotted: {dstr}'
      clrb_pos = [0.1, 0.07, 0.8, 0.04]
    else:
      btx0     = ''
      stxt     = ''
      clrb_pos = []

    mutob.plot_Nxsections(dF, Xsgm, ZM, Hbtm, Xbtm, clrmp, \
              rmin, rmax, xl1, xl2, fig1, ax1, sttl=sttl, stxt=stxt,\
              btx=btx0, clrb_ornt='horiz', \
              clrb_pos=[0.1, 0.07, 0.8, 0.02], txt_pos=[0.1, 0.9, 0.8, 0.09])



