"""
  Check relax fields from PIOMAS monthly ice thickness and concentration
  see: piomas_relaxation_yearly.py

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

YR0 = 1993
MM0 = 4
ifld = 'ithkn'  # ithkn, iarea
file_type = 'monthly'  # monthly, daily, ... or clim
                       # for climatologies, do not need padded time - data will be recycled
                       # for monthly, daily, etc. need -dt and +dt at the beginn/end 
plot_test_tile = True  # show tile with the test grid point that print statistics during SIS2 run

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
flout = f'PIOMAS_ithkn_iconc_{YR0}_{file_type}.nc'
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

AA = ds_rlx[ifld].isel(time=itime).data
AA = np.where(HH>=0, np.nan, AA)

# Global indices of the tile:
isdG = 305
iedG = 323
jsdG = 674
jedG = 691
itestG = 317
jtestG = 684
i0 = itestG-1
j0 = jtestG-1
print(f"test pnt i/j = {i0}/{j0} {ifld}={AA[j0,i0]:.3f}")

# Read saved grid indices and ice fields for test tile:
pthtxt = '/work/Dmitry.Dukhovskoy/run_output/NEP_ISPONGE/1993/04'
ftxt = os.path.join(pthtxt,'test_sponge.txt')

II = []
JJ = []
VALS = []
with open(ftxt, 'r') as file:
  for line in file:
    # Process each line
    astr= line.strip()
    strl = astr.split('  ')

    ii = int(strl[0]) - 1 # convert to python
    jj = int(strl[1]) - 1 # convert to python
    val = float(strl[2])

    II.append(ii)
    JJ.append(jj)
    VALS.append(val)

II = np.array(II)
JJ = np.array(JJ)
VALS = np.array(VALS)

# Plot all points with 0 values:
indx0 = np.where(VALS <= 1.e-10)[0]

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


sttl = f'Relaxation {ifld} SIS2 from PIOMAS {YR0}/{MM0}'

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
img = ax1.pcolormesh(AA, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.plot([isdG,iedG],[jsdG,jsdG],'-')
ax1.plot([isdG,iedG],[jedG,jedG],'-')
ax1.plot([iedG,iedG],[jsdG,jedG],'-')
ax1.plot([isdG,isdG],[jsdG,jedG],'-')

ax1.plot(II[indx0],JJ[indx0],'.')


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

btx = 'check_test_tile.py' 
bottom_text(btx, pos=[0.2, 0.01])


