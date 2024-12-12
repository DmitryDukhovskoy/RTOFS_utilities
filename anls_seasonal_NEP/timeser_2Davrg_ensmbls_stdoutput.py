"""
  Plot variables averaged over some area 
  for ensemble runs

  Use standard output daily files
  ocean_daily.nc

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
importlib.reload(mutob)

# experiment: year start, month start, ...
# change dayrun to plot desired date output - # of days since start date
# variables in ocean_daily.nc:
# sos - Sea Surface Salinity
# ssh 
# tos - Sea Surface Temperature
# tob - Sea Water Potential Temperature at Sea Floor
# sob - Sea Water Salinity at Sea Floor
# For bottom variables - look at shelf region only <400m depth
#
varnm  = 'sob'  # temp (potential) / salin
hmin   = -400.

expt_nmb = 2    # =1 - OBs from fixed SPEAR ens #, =2 - OBs from multi-ens. SPEAR 
YRS    = 1993 # year start of the forecast
MOS    = 4
DDS    = 1    
nensR  = 1  #ens # for reference ensemble run
regn   = 'poly_north'  # region to do the averaging over
archv_fl = 'ocean_daily.nc' # output file name 
lr     = -1             # vertical layer to analyze, =-1 for 2D output fields

expt     = "seasonal_daily"  # seas f/casts with dailyOB from SPEAR
runname  = f'NEPphys_frcst_dailyOB-expt{expt_nmb:02d}'
#'_{YRS}-{MOS:02d}-e{nensR:02d}'
dnmbS    = mtime.datenum([YRS,MOS,DDS]) 
dv_start = mtime.datevec(dnmbS)

print(f'Expt: {expt} Run: {runname} init date: {YRS}/{MOS}/{DDS}')

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

pthoutp    = pthseas['MOM6_NEP'][expt]['pthoutp'].format(expt_nmb=expt_nmb)
#pthwoutp   = pthseas['MOM6_NEP'][expt]['pthwoutp'].format(runname=runname)
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

# Select some region:
# Get indices of the polygon:
II = pthseas['ANLS_NEP'][regn]['II']
JJ = pthseas['ANLS_NEP'][regn]['JJ']
jdm, idm = HH.shape

DX, DY = mmom6.dx_dy(hlon, hlat)
Acell  = DX*DY
X, Y   = np.meshgrid(np.arange(idm), np.arange(jdm))
MS, _, _ = mmisc.inpolygon_v2(X, Y, II, JJ)  # 
if varnm == 'tob' or varnm == 'sob':
  JBS, IBS = np.where( (MS == 1) & (HH < 0) & (HH > hmin) ) #exclude deeep regions
else:
  JBS, IBS = np.where( (MS == 1) & (HH < 0) ) #exclude deeep regions
MSKBS  = np.zeros((jdm,idm))
MSKBS[JBS,IBS] = 1

ENSR = [x for x in range(1,11)]
NensR = len(ENSR)
for iens in range(NensR):
  nens = ENSR[iens]
  pthfcst   = os.path.join(pthoutp,f"{YRS}-{MOS:02d}-e{nens:02d}","history")

  print(f"Processing ens={nens:02d} {pthfcst}")

  Fts, TM = manseas.timeser_spatavrg_stdoutp(pthfcst, archv_fl, varnm, lr, MSKBS, Acell)

#  dF = F2d - F2dR
#  print(f"diff min/max: {np.nanmin(dF)}/{np.nanmax(dF)}")
  if iens == 0:
    dim1 = "Time"
    darr_var = xarray.DataArray(Fts, dims=(dim1), coords={dim1: TM})
    dset1D = xarray.Dataset({f"{varnm}_e{nens:02d}": darr_var}).copy()
  else:
    darr_var = xarray.DataArray(Fts, dims=(dim1), coords={dim1: TM})
    dset_var = xarray.Dataset({f"{varnm}_e{nens:02d}": darr_var})
    dset1D = xarray.merge([dset1D, dset_var])

# 
match regn:
  case "poly_south":
    if varnm == 'tos':
      yl1, yl2 = 17.5, 27.0
    elif varnm == 'sos':
      yl1, yl2 = 33.95, 34.3
    elif varnm == 'ssh':
      yl1, yl2 = 0.1, 0.6
    elif varnm == 'tob':
      yl1, yl2 = 12.5, 17.8
    elif varnm == 'sob':
      yl1, yl2 = 34.2, 34.45

  case "poly_central":
    if varnm == 'tos':
      yl1, yl2 = 17.5, 27.0
    elif varnm == 'sos':
      yl1, yl2 = 32.0, 32.8
    elif varnm == 'ssh':
      yl1, yl2 = 0., 0.35
    elif varnm == 'tob':
      yl1, yl2 = 7.5, 11.
    elif varnm == 'sob':
      yl1, yl2 = 32., 34.

  case "poly_north":
    if varnm == 'tos':
      yl1, yl2 = 10.5, 27.0
    elif varnm == 'sos':
      yl1, yl2 = 31.75, 32.6
    elif varnm == 'ssh':
      yl1, yl2 = -0.15, 0.15
    elif varnm == 'tob':
#      yl1, yl2 = 1.6, 2.0
      yl1, yl2 = 5.0, 8.5
    elif varnm == 'sob':
#      yl1, yl2 = 34.4, 34.52 
      yl1, yl2 = 32.3, 33.0

# ===================
# Plotting
# ===================
#btx = 'timeser_2Davrg_ensmbl.py' 
btx = 'timeser_2Davrg_ensmbls_stdoutput.py' 

plt.ion()

fgnmb = 1
fig1 = plt.figure(fgnmb,figsize=(9,8))
plt.clf()

ax1 = plt.axes([0.1, 0.34, 0.8, 0.6])
LNS  = []
for ie in range(NensR):
  nens = ENSR[ie]
  vards = f"{varnm}_e{nens:02d}"
  Ts = dset1D[vards].data.squeeze()
  TM = dset1D['Time'].data
  Tday = TM-TM[0]

  ln1, = ax1.plot(Tday, Ts, '-', label=f'ens{nens:02d}')
  LNS.append(ln1)


ax1.set_ylim([yl1,yl2])
ax1.grid('on')
ax1.set_xlabel('Forecast days')
ax1.set_ylabel(f'{varnm}')


if varnm == 'salin': ax1.set_ylabel('S')
dstart = f'{dv_start[0]}/{dv_start[1]}/{dv_start[2]}'
sttl = f'{runname} for bottom T,S: hmin={hmin:.0f}m\n Seas f/cast init: {dstart}, {varnm} {archv_fl} {regn}'
ax1.set_title(sttl)


ax3 = plt.axes([0.55, 0.04, 0.4, 0.27])
lgd = plt.legend(handles=LNS, loc='upper right')
ax3.axis('off')


bottom_text(btx)


f_showreg = False
if f_showreg:
  dfmom6 = os.path.join(pthfcst, archv_fl)
  dset   = xarray.open_dataset(dfmom6)

  A2d = dset[varnm].isel(time=30).data.squeeze()

# Stereographic Map projection:
  from mpl_toolkits.basemap import Basemap, cm
  lon0 = 220.
  lat0 = 60.
  res  = 'l'
  m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution=res)
  xR, yR = m(hlon, hlat)
  PMsk = ( (xR > 1e20) | (yR > 1e20) )
  AA = A2d.copy()
  AA = np.insert(AA, 0, AA[:,-1], axis=1)
  AA = np.insert(AA, -1, AA[-1,:], axis=0)
  AA[PMsk] = np.nan
  xR[PMsk]   = 1.e30
  yR[PMsk]   = 1.e30

  ny, nx = HHG.shape
  AA = AA[0:ny, 0:nx]
  AA = np.where(HHG >= 0., np.nan, AA)


  xBND, yBND = m(II, JJ)




  m = Basemap(width=3300*1.e3,height=3300*1.e3, resolution='l',\
              projection='stere', lat_ts=55, lat_0=62, lon_0=-175)

  xR, yR = m(hlon, hlat)

  xcMap, ycMap = m(hlon[JJ,II], hlat[JJ,II])

  plt.ion()

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  m.drawcoastlines()
  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

  img = m.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  m.plot(xcMap, ycMap, '-', color=[0.8, 0.4, 0])
#  m.contour(xR, yR, HH, [-1000], colors=[(0,0,0)], linestyles='solid')
#  m.contour(xR, yR, Tbtm, [2.], colors=[(0,0.5,1.)], linestyles='solid')
  #ax1.axis('scaled')
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

  ax3 = fig1.add_axes([0.02, 0.02, 0.8, 0.05])
  ax3.text(0, 0, sinfo, fontsize=10)
  ax3.axis('off')

  btx = 'plot_coldpool.py'
  bottom_text(btx, pos=[0.2, 0.01])


# Bind the button_press_event with the onclick() method
# Function to print mouse click event coordinates
  def onclick(event):
     print([event.xdata, event.ydata])

  fig1.canvas.mpl_connect('button_press_event', onclick)

  fig1.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  ax1.contour(HH,[0], linestyles='solid', colors=[(0., 0.2, 0.5)])
  ax1.contour(HH,[-1000], linestyles='solid', colors=[(0., 0.6, 0.9)])
  ax1.axis('scaled')




