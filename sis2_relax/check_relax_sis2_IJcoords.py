"""
  Check relax fields from PIOMAS monthly ice thickness and concentration
  Plot in I-J coords for easy checking

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
import pickle
from yaml import safe_load

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

plot_fields = True
plot_piomas = False
YR1 = 1993
YR2 = 1994 
YR0 = 1993   # year to plot
MM0 = 6      # month to plot
ifld = 'iarea'  # ithkn, iarea
# Test point in Fortran indices:
iF0 = 230  
jF0 = 700
i0 = iF0-1 ; j0 = jF0-1
file_type = 'monthly'  # monthly, daily, ... or clim
                       # for climatologies, do not need padded time - data will be recycled
                       # for monthly, daily, etc. need -dt and +dt at the beginn/end 

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

# Read saved relax. fields:
flout = f'PIOMAS_ithkn_iconc_{YR1}_{YR2}_{file_type}.nc'
diclim = os.path.join(pthsis, flout)
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

# Read PIOMAS field:
match ifld:
  case('ithkn'):
    varnm = varthck
    dfpiomas = os.path.join(pthdata,flthck)
    clrmp = mclrmps.colormap_ice_thkn()
    rmin = 0.
    rmax = 4.
  case('iarea'):
    varnm = varconc
    dfpiomas = os.path.join(pthdata,flconc)
    clrmp = mclrmps.colormap_conc()
    rmin = 0.
    rmax = 1.

clrmp.set_bad(color=[0.2, 0.2, 0.2])
A2dP = msisrlx.read_PIOMAS(YR0, MM0, dfpiomas, varnm)

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])

sttl = f'Relaxation {ifld} SIS2 from PIOMAS {YR0}/{MM0}\n'
sttl = sttl + f"Test pnt iF0/jF0 = {iF0}/{jF0}, year={YR0}, MM0={MM0}, {ifld}: {A2dS[j0,i0]:.6f}"

img = ax1.pcolormesh(A2dS, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.contour(HH,[0], linestyles='solid', linewidths=1, colors=[(0.8, 0.8, 0.8)])
ax1.plot(i0,j0,'o')

ax1.axis('scaled')
ax1.set_xlim([100,idm])
ax1.set_ylim([550,jdm])
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







