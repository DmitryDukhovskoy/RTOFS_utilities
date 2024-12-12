"""
  Plot 2D mean T/S fields to analyze seasonal water mass structure
  in different regions: # CalCur - Calif Current region, Alaska - Alaska region, BerSea - Bering
  Following Stoke et al., 2015

  in the Calif. Current region, see analysis:
  Auad et al., 2011
  The California Current System in relation to the Northeast Pacific Ocean circulation
  https://www.sciencedirect.com/science/article/pii/S0079661111001157

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
import pickle
from copy import copy
import matplotlib.colors as colors
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
sys.path.append(PPTHN + '/TEOS_10/gsw')
sys.path.append(PPTHN + '/TEOS_10/gsw/gibbs')
sys.path.append(PPTHN + '/TEOS_10/gsw/utilities')
import mod_swstate as msw
import conversions as gsw

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
importlib.reload(mutob)


# Initial date
# Look at ens run #1 - the only ens. that has 5-day av. output fields
expt     = 'seasonal_daily'  # seasonal forecasts with dailyOB from SPEAR
varnm    = 'salin'  # temp (potential) / salin
#dnmbS    = mtime.datenum([2015,1,1])
# Averaging time period:
MMS   = 1    # init. month in each year, can be changed to different months: 1, 4, 7, 10
YAVRG = [x for x in range(2011,2021)]
MAVRG = [1,2,3]  # months to average:
regn_name = 'CalCur' # CalCur - Calif Current region, Alaska - Alaska region, BerSea - Bering
                     # Following Stoke et al., 2015

nensR    = 1
expt_nmb = 2   # 2 - seas f/casts with dailyOB


lr0  = 1  # ocean layers from 1, ..., 50

expt_name = f'NEPphys_frcst_climOB{expt_nmb:02d}'
run_info = f'{expt_name} init MM={MMS} e{nensR:02d}, avrg: Years:{min(YAVRG)}-{max(YAVRG)} Mo: {min(MAVRG)}-{max(MAVRG)}'

print(f'Plotting {varnm} {expt_name} ')
print(f'{run_info}')

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
ndav       = pthseas['MOM6_NEP'][expt]['ndav']  # # of days output averaged

# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

#mo_fcsts = manseas.yrmo_seasonal_fcst(YRS, MMS)

ocnfld = 'oceanm'
pthoutp0 = pthseas['MOM6_NEP'][expt]['pthoutp'].format(expt_nmb=expt_nmb)

Time = []
iyr  = 0
for YRS in (YAVRG):
  dnmbS    = mtime.datenum([YRS,MMS,1])
  pthfcst0 = os.path.join(pthoutp0,f'{YRS}-{MMS:02d}-e{nensR:02d}','history')
  AA, TM = manseas.monthly_mean_from_Ndaily_ocean2D(pthfcst0, YRS, MMS, varnm, ocnfld, lr0, MAVRG=MAVRG)

  if iyr == 0:
    A2d = AA.copy()
  else:
    A2d = A2d + AA
  Time = Time + TM

  iyr += 1

A2d = A2d/iyr
DV = mtime.datevec2D(Time)

II = pthseas['ANLS_NEP'][regn_name]['II']
JJ = pthseas['ANLS_NEP'][regn_name]['JJ']
xlim1 = min(II)
xlim2 = max(II)
ylim1 = min(JJ)
ylim2 = max(JJ)


if varnm == 'salin' or varnm == 'salt': 
  clrmp = mclrmps.colormap_haline2()
  clrmp.set_bad(color=[0., 0., 0.])
  rmin = 30.0
  rmax = 35.0
  tscntrs = [x/10 for x in range(320,360,2)]
  tslabels = [x for x in range(32,36)]
  cntr_clr = [0.3, 0.3, 0.3]

elif varnm == 'temp' or varnm == 'potT': 
  clrmp = mutil.colormap_temp(clr_ramp=[0.9,0.8,1])
  clrmp.set_bad(color=[1,1,1])
  rmin = -2.
  rmax = 23.
elif varnm == 'ssh':
  clrmp = mutil.colormap_ssh(nclrs=200)
  rmin = -0.5
  rmax = 0.5


plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
im1 = ax1.pcolormesh(A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
ax1.axis('scaled')
ax1.set_xlim([xlim1, xlim2])
ax1.set_ylim([ylim1, ylim2])
ax1.contour(hlon, [x for x in range(200, 360, 10)], colors=[(0.8, 0.8, 0.8)], linestyles='solid', linewidths=1)
ax1.contour(hlat, [x for x in range(0, 89, 10)], colors=[(0.8, 0.8, 0.8)], linestyles='solid', linewidths=1)

CS = ax1.contour(A2d, tscntrs, colors=[cntr_clr], linestyles='solid', linewidths=1)
ax1.clabel(CS, tslabels, inline=1, fontsize=10)

sttl = f"{run_info}"
ax1.set_title(sttl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
# extend: min, max, both
clb = plt.colorbar(im1, cax=ax2, orientation='vertical', extend='both')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)


btx = 'plot_timeavrgTS_regions.py'
bottom_text(btx, fsz=6, pos=[0.05, 0.03])














