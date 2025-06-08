"""
  Some ERA5 padded forcing files cause ice continuity 
  blow up (negative ice thickness)

  Happened on Dec. 31 2000
  
  Try to modify the last day (padding) of the ERA5 padded files
  Replace the day following the model blow up by copying
  the fields at -dltHR over to N hours/records of data into all hrly slots


#   1-hr fields are assumed
#
# ---*----X---*---*---X---*---*---*---*
#   time-dltHR      date       time+dltHR
#                 to replace
#
#         |---------->|   Data from (time-dltHR) copied
#

  ERA5 atm. forcing files:
  ERA5_u10_2000_padded.nc
  ERA5_v10_2000_padded.nc
  ERA5_lp_2000_padded.nc
  ERA5_msl_2000_padded.nc
  ERA5_sf_2000_padded.nc
  ERA5_sphum_2000_padded.nc
  ERA5_ssrd_2000_padded.nc
  ERA5_strd_2000_padded.nc
  ERA5_t2m_2000_padded.nc
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
import pandas as pd

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
import mod_interp1D as mint1d
#importlib.reload(mutob)

parser = argparse.ArgumentParser()
parser.add_argument("--yr", help="year of file to change: 1993, ..., 2020", type=int)
parser.add_argument("--date", help="Date needs to be changed: YYYYMMDD", type=int)
parser.add_argument("--dlthr",help="Time window (hours) +/- from date to do time interpolation", type=int)
parser.add_argument("--varnm", help="field to change: u10,v10,lp,msl,sf,sphum,ssrd,strd,t2m", type=str)
args = parser.parse_args()

#VARS = ['u10','v10','lp','msl','sf','sphum','ssrd','strd','t2m']
#VARS = ['msl','sf','sphum']
#VARS = ['ssrd','strd','t2m']
VARS = ['u10','v10','lp']

f_save = False

dltHR = 1.
if args.varnm:
  VARS = [args.varnm]
if args.yr:
  YR = args.yr
if args.date:
  date_change = int(args.date)
  yr_ch, mm_ch, dd_ch = mtime.extract_yymmdd(date_change)
  dnmb_ch = mtime.datenum([yr_ch,mm_ch,dd_ch])
if args.dlthr:
  dltHR = args.dlthr

# Start-end dates to replace field:
dnmb_prev = dnmb_ch-dltHR/24.0
dnmb_next = dnmb_ch+dltHR/24.0

if not f_save:
  print(f'Files not saved at the end, f_save flag is off ...')


pthoutp = f'/work/Dmitry.Dukhovskoy/NEP_input/ERA5_padded_changed/{YR}'
os.makedirs(pthoutp, exist_ok=True)

def read_ERA5_field(varnm,dnmbR):
  """
    Read ERA5 field for given date dnmbR
    Assumed hourly data
  """
  YRR,MMR,DDR,hrr,minr = mtime.datevec(dnmbR, round_hrs=True)

  pthera = '/archive/e1n/mom6/NEP/atmos_forcing/era5_padded'
  ptherafld = os.path.join(pthera,varnm)
  flera = f'ERA5_{varnm}_{YRR}_padded.nc'
  dflera = os.path.join(ptherafld,flera)

  print(f'Processing {YRR} {varnm}')
  print(f'Opening {dflera}')
  # ds = xarray.open_dataset(dflera, mask_and_scale=False)  # keep missing values unmasked
  dset = xarray.open_dataset(dflera)
  Time = dset['time'].data
  tmP = pd.to_datetime(Time)
  nrec = len(tmP)
  TNEP = np.zeros((nrec,4), dtype=int)
  years  = tmP.year.to_numpy()
  months = tmP.month.to_numpy()
  days   = tmP.day.to_numpy()
  hours  = tmP.hour.to_numpy()
  TM = np.zeros((nrec))
  #TM     = mtime.datenum([years,months,days,hours]) <-- need to change mtime.datenum to work with 1D arrays
  for irec in range(nrec):
    yy,mm,dd,hh = years[irec],months[irec],days[irec],hours[irec]
    TM[irec] = mtime.datenum([yy,mm,dd,hh])

  # Check if requested date is in the time range:
  assert(dnmbR >= TM[0] and dnmbR <= TM[-1]), 'Requested date is outside the time range in the file'    
  DD = abs(TM-dnmbR)
  itime = np.argmin(DD)
  assert(DD[itime] < 1.e-3),f'Could not find requested date  {dnmbR}'
  A2d = dset[varnm].isel(time=itime).data

  return A2d, TM, dset
  
ich = 181
jch = 237
nvars = len(VARS)
for varnm in VARS:
  # Read previous date from the data set being modified:
  Aprev, TMp, dset = read_ERA5_field(varnm, dnmb_prev)
  darray = dset[varnm]
  dmm1 = darray.data[:,jch,ich].copy()

  # Define Time array of time being replaced:
  # Hourly data assumed
  # Cannot go beyond the last padded date
  dltT = 1.0 / 24.0  # 1 hour in days
  Tend = np.min([dnmb_next - dltT, TMp[-1]])
  Tint = np.arange(dnmb_prev + dltT, Tend + dltT / 2, dltT)
  yrP,mmP,ddP,hrP = mtime.datevec(dnmb_prev, round_hrs=True)[:4]

  jdm,idm = Aprev.shape
  for it in range(len(Tint)):
    dnmbI = Tint[it]
    yri,mmi,ddi,hri = mtime.datevec(dnmbI, round_hrs=True)[:4]
    print(f'Copying data {yrP}/{mmP}/{ddP}:{hrP:02d} --> {yri}/{mmi}/{ddi}:{hri:02d}')  

    # Replace with the reference field:
    DT = abs(TMp-dnmbI)
    itime = np.argmin(DT)
    assert(DT[itime]<1./24.),f'Time differes > 1hr, {dnmbI:.5f} vs {TMp[itime]:.5f}' 
    # CHeck introduced difference:
    D2 = darray.data[itime,:,:]-Aprev
    rmax = np.max(D2)
    rmin = np.min(D2)
    print(f"min/max error after interpolation = {rmin:.4f}/{rmax:.4f}")
    Dold = darray.data[itime,:,:].copy()    # for debugging
    darray.data[itime,:,:] = Aprev

  dmm2 = darray.data[:,jch,ich].copy()
  dset[varnm] = darray
  dset.attrs["info"]=f"Modified {yr_ch}/{mm_ch}/{dd_ch}:0hr - {dltHR}hrs "
  dset.attrs["code"]="/home/Dmitry.Dukhovskoy/python/sis2_relax/era5_change_copyNrecs.py"

  if f_save:
    flnew = f'ERA5_{varnm}_{YR}_cp{int(dltHR):03d}hrs_padded.nc'
    dflout = os.path.join(pthoutp, flnew)
    print(f"Saving to {dflout}")
    dset.to_netcdf(
      dflout,
      format='NETCDF3_64BIT',
      engine='netcdf4',
    )     

f_chck = False
if f_chck:
  clrmp_dlt = mclrmps.colormap_ssh(cpos='YlOrRd', cneg='PuBu_r')
  clrmp_dlt.set_bad(color=[0.2,0.2,0.2])
  dmin = -5.
  dmax = 5.

  dI = Dold-Aprev
  
  yrP,mmP,ddP,hrP = mtime.datevec(dnmb_prev, round_hrs=True)[:4]
  sttl = f'Difference {varnm} {yri}/{mmi}/{ddi}:{hri:02d}h - {yrP}/{mmP}/{ddP}:{hrP:02d}h'

  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  img = ax1.pcolormesh(dI, cmap=clrmp_dlt, vmin=dmin, vmax=dmax)
  ax1.set_title(sttl)

  ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                     0.02, ax1.get_position().height])
  # extend: min, max, both
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
  ax2.yaxis.set_ticks(list(np.linspace(dmin,dmax,11)))
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)









