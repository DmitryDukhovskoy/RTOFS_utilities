"""
  Interpolate NSIDC sea ice concentration 
  to NEP grid
  gmapi indices: get_gmapi_NSIDC_to_SIS2.py

  Need to update older version 3:
  Latest version of NSDIC NRT v5 is available here:
  https://noaadata.apps.nsidc.org/NOAA/G02202_V5/north/monthly/
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
import mod_plot_xsections as mxsct
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
#import mod_valid_utils as mvutil
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

parser = argparse.ArgumentParser()
parser.add_argument("--YRS", help="year to start interp.: 1993, ..., 2020", type=int)
parser.add_argument("--YRE", help="year to end interp.: 1993, ..., 2020", type=int)
#parser.add_argument("--MMS", help="month to start interp, defualt=1 : 1,..., 12", type=int)
#parser.add_argument("--MME", help="month to end interp, default=12 or =MMS: 1,..., 12", type=int)
args = parser.parse_args()

# experiment: year start, month start, ...
# change dayrun to plot desired date output - # of days since start date
# in daily-mean output fields: date is in the middle of the averaging period
varnm  = 'iconc'

# Default values that can be modified by keywords
# Day to plot either in year days or actual date:
YRS  = 1993
YRE  = YRS
MMS  = 1
MME  = 12
f_save = False

if args.YRS:
  YRS = args.YRS
  YRE = YRS
if args.YRE:
  YRE = args.YRE


fyaml_param='relax_expts.yaml'
with open(fyaml_param) as ff:
  param_expt = safe_load(ff)  

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

expt     = "seasonal_daily"
pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape
LMsk = np.where(HH<0, 1, 0)

# Mask out southern lats:
LMsk[:551,:] = 0
LMsk[:,:47] = 0


# Get gmapi 4 NSIDC grid points for interpolation
pthdump = pthseas['NRT_NSIDC']['pthdump']
fgmapi  = f'NSIDC_NRTice_NEP_gmapi_{jdm}x{idm}.nc'
dfgmapi = os.path.join(pthdump, fgmapi)
print(f'Loading gmapi --> {dfgmapi}')

dgmapi = xarray.open_dataset(dfgmapi)
IMOM = dgmapi['mom_indx'].data
JMOM = dgmapi['mom_jndx'].data
INDX = dgmapi['gmapi_i'].data
JNDX = dgmapi['gmapi_j'].data

icc = 0
LON = None
LAT = None
mcal = np.arange(1,13)
for YR in range(YRS,YRE+1):
  A3d = np.zeros((12,jdm,idm))
  for MM in mcal:
    print(f'Processing {YR}/{MM}')
    # Note that data fields are flipped upside-down wrt the coord. grid
    # shown here:
    # https://nsidc.org/data/user-resources/help-center/guide-nsidcs-polar-stereographic-projection
    # this does not correspond to the X/Y coordinates written in the netcdf files
    # see: Xnrt  = dset_nsidc['xgrid'].data
    # Ynrt  = dset_nsidc['ygrid'].data
    #
    # gmapi are computed for the grid using X/Y metric coordinates
    # Need to flip the data to use these gmapi
    fsfx = 'f11'
    if (YR == 1995 and MM >= 10) or (YR > 1995 and YR <2008):
      fsfx = 'f13'
    if (YR >= 2008):
      fsfx = 'f17'

    #pthnsidc = f'/work/Dmitry.Dukhovskoy/data/NRT_NOAA_NSIDC_seaconc/{YR}_mnth'
    pthnsidc = pthseas['ALL']['dirnsidc_intrp'].format(YR=YR)
    flnsidc = f'seaice_conc_monthly_nh_{YR}{MM:02d}_{fsfx}_v04r00.nc'
    dset = xarray.open_dataset(os.path.join(pthnsidc,flnsidc))

    if LON is None or LAT is None:
      # Reconstructed lon/lat from metric Polar Sterographic proj.
      # Note orientation of the LON/LAT vs data - do not match!
      _, LON, LAT = manseas.avrg_cice_NSIDC(2000, 2000, 1, 1)

    #
    # Choose NASA algorithm sea ice conc:
    AA = dset['nsidc_nt_seaice_conc_monthly'].data.squeeze()
    CIce = np.flipud(AA)
    CIce = np.where(CIce>2., np.nan, CIce)
    CIint = msisrlx.interp2Dfld(CIce, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)
    CIint = np.where(HH>=0, np.nan, CIint)
    A3d[MM-1,:,:] = CIint
  
  darr_cice = xarray.DataArray(A3d, dims=("time","jdim","idim"),\
                     coords={"time": mcal,\
                             "jdim": np.arange(jdm),\
                             "idim": np.arange(idm)})
  dset = xarray.Dataset({"ice_conc": darr_cice})
  dset['ice_conc'].attrs['long_name']='ice partial area'
  # Add global attributes:
  dset.attrs.update({
    "info": "Interpolated from NSIDC NRT ice conc. NASA retrieval algorithm",
    "code": "interp_NSIDC_NEP_mnthly.py"
  })

  if f_save:
    fliceout = f'NSIDC_iconc_mnth_interpNEP{jdm}x{idm}_{YR}.nc'
    dfliceout = os.path.join(pthnsidc,fliceout)
    print(f'Dumping interpolated ice conc --> {dfliceout}')
    dset.to_netcdf(dfliceout, format='NETCDF3_64BIT', engine='netcdf4')
  else:
    print("Not saved, turn f_save flag on")


# -------------------
#
# Plot ice fields
#
# -------------------
f_plot = False
if f_plot:
  MM = 6
  imo = MM-1
  A2d = dset['ice_conc'].isel(time=imo).data
 
  clrmp = mclrmps.colormap_conc()
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 0.
  rmax = 1.
  if varnm == 'ithkn':
    clrmp = mclrmps.colormap_ice_thkn()
    clrmp.set_bad(color=[0.2, 0.2, 0.2])
    rmin = 0.
    rmax = 4.

  sttl = f"NSIDC NRT ice conc, {YR}/{MM:02d}"
  # Stereographic Map projection:
  from mpl_toolkits.basemap import Basemap, cm
  m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
              projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

  xR, yR = m(hlon, hlat)


  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  m.drawcoastlines()
  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

  img = m.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  #m.contour(xR, yR, HH, [-1000], colors=[(0,0,0)], linestyles='solid')
  #ax1.axis('scaled')
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

  sinfo='NSIDC NRT sea ice conc, NASA retrieval algorithm'
  ax3 = fig1.add_axes([0.02, 0.025, 0.8, 0.05])
  ax3.text(0, 0, sinfo, fontsize=8)
  ax3.axis('off')


  btx = 'interp_NSIDC_NEP_mnthly.py'
  bottom_text(btx, pos=[0.2, 0.01])


