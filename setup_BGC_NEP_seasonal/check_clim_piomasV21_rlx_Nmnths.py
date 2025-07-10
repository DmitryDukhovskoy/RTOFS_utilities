"""
  Plot irlx fields created by create_clim_piomasV21_rlx_yearly.py
  N months in 1 fig
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
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
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_mom6 as mmom6


parser = argparse.ArgumentParser()
parser.add_argument("--yrp", help="Year to plot", type=int, required=True)
parser.add_argument("--yrf", help="Start Year in the file name, default=YR", type=int)
parser.add_argument("--mms", help="Cal. month start to plot, default 1", type=int)
parser.add_argument("--mme", help="Cal. month end to plot, default 12", type=int)
parser.add_argument("--varnm", help="variable name to plot: iconc or ithkn", required=True, type=str)
parser.add_argument("--ncol", help="N of subplot columns, default - determined from N months", type=int)
parser.add_argument("--nrow", help="N of rows for subplots, default - determined from N months", type=int)
args = parser.parse_args()

YRP  = args.yrp if args.yrp else None  # year to plot
YRS  = args.yrf if args.yrf else YRP # make YRS=YRP-1 if need to check 2nd year in the file
MMS  = args.mms if args.mms else 1
MME  = args.mme if args.mme else 12
nmnths = MME-MMS+1
ncol = args.ncol if args.ncol else None
nrow = args.nrow if args.nrow else nmnths
varnm = args.varnm if args.varnm else None

if not ncol: 
  ncol = min([nmnths,4])
  nrow = nmnths // ncol
assert ncol*nrow == nmnths, f'Specified columns/rows ({ncol}/{nrow}) do not match N months {nmnths}'

Navrg = 5     # n years used for averaging PIOMAS fields
Nyrs_file = 2 # n of years in 1 file

print(f'Plotting Ice rlx climatology {varnm} {YRP} {MMS}:{MME}')

fyaml = 'bgc_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

pthtopo    = pthseas['MOM6_NEP']['hindcast']['pthgrid']
fgrid      = pthseas['MOM6_NEP']['hindcast']['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"]['hindcast']["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)

# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

#if varnm == 'iconc':
match varnm:
  case "iconc":
    varnc = 'iarea'
    clrmp = mclrmps.colormap_conc()
    clrmp.set_bad(color=[0.2, 0.2, 0.2])
    rmin = 0.
    rmax = 1.

    #ci0=0.15
    hcntrs = [0.15]
    cntr_clr = [1, 0.4, 0]

  case "ithkn":
    varnc = 'ithkn'
    clrmp = mclrmps.colormap_ice_thkn()
    clrmp.set_bad(color=[0.2, 0.2, 0.2])
    rmin = 0.
    rmax = 4.

    hcntrs = [x for x in range(1,10)]
    cntr_clr = [.95, 0.95, 0.95]

def plot_field(ax1, fig1, m, xR, yR, A2d, clrmp, rmin, rmax, plt_clrb, sttl=[]):
  fig1.sca(ax1)
  m.drawcoastlines()
  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

  img = m.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(sttl)

  # extend: min, max, both
  if plt_clrb:
    ax2 = fig1.add_axes([0.9,0.1,0.013,0.8])
    clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')

    ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
    ax2.set_yticklabels(ax2.get_yticks())
    ticklabs = clb.ax.get_yticklabels()
    #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
    clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
    clb.ax.tick_params(direction='in', length=12)

  return ax1

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

xR, yR = m(hlon, hlat)

plt.ion()
#plt.close('all')
#fig1, axes = plt.subplots(nrow, ncol, figsize=(15, 10))  # 3 rows, 4 columns
fig1 = plt.figure(1,figsize=(15, 10))  
fig1.clf()  # Clear the figure
axes = fig1.subplots(nrows=nrow, ncols=ncol) 

# Shift subplots to the left
# and up for colorbar and text
# also keep subplots close to each other : wspace, hspace
fig1.subplots_adjust(
    left=0.05,
    right=0.85,  # More right-side room
    top=0.95,
    bottom=0.1,
    wspace=0.05,
    hspace=0.1
)

pthrlx = '/work/Dmitry.Dukhovskoy/NEP_input/SIS2_relax'
YRE = YRS+Nyrs_file-1
flnm = f'PIOMASv21_ithkn_iconc_{YRS}_{YRE}_avrg{Navrg}yr.nc'

dcice = os.path.join(pthrlx,flnm)
print(f'Reading target relax fields: {dcice}')
ds_rlx = xarray.open_dataset(dcice)
Time = ds_rlx['time'].data
TM = mmisc.convert_nptime_to_datenum(Time)

sinfo = f'Rlx fields from PIOMAS {Navrg}-yr averaging\n' + dcice

iplt = 0
for MM in range(MMS,MME+1):
  dnmb0 = mtime.datenum([YRP,MM,15,12])
  D = abs(TM-dnmb0)
  itime = np.argmin(D)
  dv0 = mtime.datevec(TM[itime])
  assert dv0[0]==YRP, f'Requested YR={YRP}, year in rlx file={dv0[0]}'
  assert dv0[1]==MM, f'Requested month={MM}, month in rlx file={dv0[1]}'

  A2d = ds_rlx[varnc].isel(time=itime).data.squeeze()
  A2d = np.where(HH>=0, np.nan, A2d)

  irow = iplt // ncol
  icol = iplt % ncol
  iplt += 1

  sttl = f'{Navrg}yr avrg {varnm} {YRP}/{MM:02d}'

  ax1 = axes[irow, icol]
  if iplt == 1:
    plt_clrb = True
  else:
    plt_clrb = False

  ax1 = plot_field(ax1, fig1, m, xR, yR, A2d, clrmp,rmin,rmax,plt_clrb,sttl=sttl)
  CS1 = ax1.contour(xR,yR,A2d,hcntrs,linestyles='solid', colors=[cntr_clr], linewidths=1.2)
  if varnm == 'ithkn':
    ax1.clabel(CS1, inline=1, fontsize=10)

ax3 = plt.axes([0.05, 0.065, 0.6,0.03])
ax3.text(0, 0, sinfo, fontsize=8)
ax3.axis('off')
btx = 'check_clim_piomasV21_rlx_Nmnths.py'
bottom_text(btx)


