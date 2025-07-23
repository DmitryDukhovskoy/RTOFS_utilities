"""
  Check relax fields from PIOMAS monthly ice thickness and concentration
  Plot in I-J coords for easy checking

  ARC12 domain

  monthly fields
  1901 - 2010
  https://psc.apl.uw.edu/research/projects/piomas-20c/

"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
import matplotlib.pyplot as plt
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
parser.add_argument("--yrs",    help="year start to extract PIOMAS: 1993, ..., 2020", type=int, required=True)
parser.add_argument("--yre",    help="year end to extract PIOMAS: 1993, ..., 2020", type=int)
parser.add_argument("--yrplot", help="year to plot should be yrs <= yrplot <= yre, default = yrs", type=int)
parser.add_argument("--mm",     help="month to plot should be yrs <= yrplot <= yre, default = yrs", type=int)
parser.add_argument("--ifld",   help="relax field to plot: ithkn, iarea", type=str, required=True)
parser.add_argument("--nyrs",   help="number of years grouped in 1 relax. file, default=2", type=int)
args = parser.parse_args()

file_type = 'monthly'  # monthly, daily, ... or clim
                       # for climatologies, do not need padded time - data will be recycled
                       # for monthly, daily, etc. need -dt and +dt at the beginn/end 

YR1  = args.yrs if args.yrs else None
nyrs = args.nyrs if args.nyrs else 2
YR2  = args.yre if args.yre else YR1+nyrs-1
YR0  = args.yrplot if args.yrplot else YR1
MM0  = args.mm if args.mm else None
ifld = args.ifld if args.ifld else None

plot_fields = True

#ifld = 'iarea'  # ithkn, iarea
# Test point in Fortran indices:
iF0 = 280 
jF0 = 645
i0 = iF0-1 ; j0 = jF0-1
file_type = 'monthly'  # monthly, daily, ... or clim
                       # for climatologies, do not need padded time - data will be recycled
                       # for monthly, daily, etc. need -dt and +dt at the beginn/end 

# ARC12 grid:
ptharc  = '/work/Dmitry.Dukhovskoy/ARC12/topo_grid'
dflarc  = os.path.join(ptharc,'ocean_hgrid.nc')
dfltopo = os.path.join(ptharc,'ocean_topog.nc')

ds_topo = xarray.open_dataset(dfltopo)
HH = -(ds_topo['depth'].data)
jdm, idm = HH.shape

assert HH[300,200] < 0., f'Check sign of topography, ocean pnts should be < 0'

hlon, hlat = mmom6.read_mom6grid(dflarc, grdpnt='hgrid')

if ifld == 'iconc': ifld = 'iarea'

# Read saved relax. fields:
pthsis = '/work/Dmitry.Dukhovskoy/ARC12/irlx'
flout = f'PIOMASv21_ARC12_ithkn_iconc_{YR1}_{YR2}_{file_type}.nc'
diclim = os.path.join(pthsis, flout)
print(f"Reading irlx fields from {diclim}")
ds_rlx = xarray.open_dataset(diclim)
Time = ds_rlx['time'].data
TM = mmisc.convert_nptime_to_datenum(Time)
dnmb0 = mtime.datenum([YR0,MM0,15,12])
D = abs(TM-dnmb0)
itime = np.argmin(D)
dv0 = mtime.datevec(TM[itime])
assert dv0[0]==YR0, f'Requested YR={YR0}, year in rlx file={dv0[0]}'
assert dv0[1]==MM0, f'Requested month={MM0}, month in rlx file={dv0[1]}'

A2dS = ds_rlx[ifld].isel(time=itime).data
#A2dS = np.where(HH>=0, np.nan, A2dS)

print(f"Test pnt i/j = {i0}/{j0}, year={YR0}, MM0={MM0}, {ifld}: {A2dS[j0,i0]:.6f}")

match ifld:
  case('ithkn'):
    clrmp = mclrmps.colormap_ice_thkn()
    rmin = 0.
    rmax = 4.
  case('iarea'):
    clrmp = mclrmps.colormap_conc()
    rmin = 0.
    rmax = 1.

clrmp.set_bad(color=[0.2, 0.2, 0.2])


plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])

sttl = f'Relaxation {ifld} SIS2 from PIOMAS {YR0}/{MM0}\n'
sttl = sttl + f"Test pnt iF0/jF0 = {iF0}/{jF0}, year={YR0}, MM0={MM0}, {ifld}: {A2dS[j0,i0]:.6f}"

img = ax1.pcolormesh(A2dS, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.contour(HH,[0], linestyles='solid', linewidths=1, colors=[(0.8, 0.8, 0.8)])
ax1.plot(i0,j0,'o')

plt_lonlat = True
if plt_lonlat:
  clon = [x-180 for x in range(0,360,10)]
  clat = [x for x in range(40,89,10)]

  ax1.contour(hlon, clon, linestyles='solid', colors=[(0.9,0.9,0.9)])
  ax1.contour(hlat, clat, linestyles='solid', colors=[(0.9,0.9,0.9)])


ax1.axis('scaled')
#ax1.set_xlim([100,idm])
#ax1.set_ylim([550,jdm])
ax1.set_title(sttl)

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

btx = 'check_relax_sis2_IJcoords.py' 
bottom_text(btx, pos=[0.2, 0.01])







