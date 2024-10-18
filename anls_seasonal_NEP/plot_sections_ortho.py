"""
  Show sections on orthographic projection
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
import mod_mom6 as mmom6
import mod_utils_ob as mutob
importlib.reload(mutob)

# experiment: year start, month start, ...
# change dayrun to plot desired date output - # of days since start date
# in daily-mean output fields: date is in the middle of the averaging period
sctnm = 'xsct_1000m'
YRS   = 1993
MOS   = 4
nens  = 2
expt  = f'NEPphys_frcst_climOB_{YRS}-{MOS:02d}-e{nens:02d}'

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

# Global TOPO GOFS:
import mod_read_hycom as mhycom
pthgrid      = gridfls["GOFS3.1"]["93.0"]["pthgrid"]
ftopo        = gridfls["GOFS3.1"]["93.0"]["ftopo"]
fgrid        = gridfls["GOFS3.1"]["93.0"]["fgrid"]
LONG, LATG, HHG = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)


pthfcst    = pthseas['MOM6_NEP']['seasonal_fcst']['pthoutp'].format(runname=expt)
pthtopo    = gridfls['MOM6_NEP']['seasonal_fcst']['pthgrid']
fgrid      = gridfls['MOM6_NEP']['seasonal_fcst']['fgrid']
ftopo_mom  = gridfls["MOM6_NEP"]["seasonal_fcst"]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)

# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')


def get_section(sctnm):
  Is  = pthseas['ANLS_NEP'][sctnm]['II']
  Js  = pthseas['ANLS_NEP'][sctnm]['JJ']
  IJ     = np.column_stack((Is, Js))
  DX, DY = mmom6.dx_dy(hlon,hlat)
  SGMT   = mmisc.define_segments(IJ, DX, DY, curve_ornt='positive', check_pole=False)
  II     = SGMT.I_indx
  JJ     = SGMT.J_indx
  XX     = hlon[JJ,II]
  YY     = hlat[JJ,II]
  nLeg   = SGMT.Leg_number
  hLsgm1 = SGMT.half_Lsgm1
  hLsgm2 = SGMT.half_Lsgm2
  LSgm   = np.zeros((len(II)))  # total segment length = half1 + half2
  for ik in range(len(II)):
     LSgm[ik] = hLsgm1[ik] + hLsgm2[ik]

  return II, JJ, XX, YY, LSgm


# Set up orthographic projection
from mpl_toolkits.basemap import Basemap, cm

# Add extra row/col for plotting
# Add extra row/col for plotting with pcolormesh
lonw = LONG.copy()
latw = LATG.copy()
lonw = np.insert(lonw, -1, lonw[:,-1]+0.01, axis=1)
lonw = np.insert(lonw, -1, lonw[-1,:]+0.01, axis=0)
latw = np.insert(latw, -1, latw[-1,:]+0.01, axis=0)
latw = np.insert(latw, -1, latw[:,-1]+0.01, axis=1)

II, JJ, XX, YY, LSgm = get_section(sctnm)
if np.max(YY) < 50:
  lon0 = 220.
  lat0 = 40.
else:
  lon0 = 220.
  lat0 = 60.

res  = 'l'
m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution=res)
xR, yR = m(lonw,latw)
PMsk = ( (xR > 1e20) | (yR > 1e20) )
AA = HHG.copy()
AA = np.insert(AA, 0, AA[:,-1], axis=1)
AA = np.insert(AA, -1, AA[-1,:], axis=0)
AA[PMsk] = np.nan
xR[PMsk]   = 1.e30
yR[PMsk]   = 1.e30

ny, nx = HHG.shape
AA = AA[0:ny, 0:nx]
AA = np.where(HHG >= 0., np.nan, AA)

# Find section in orth coordinates:
xSCT = []
ySCT = []
xSCT, ySCT = m(XX, YY)

# Distance along the section
# normalize by the total distance
Lsection = mmisc.dist_sphcrd(YY[-1],XX[-1],YY[0],XX[0]) # total length of the section, m
Xdist = np.cumsum(LSgm)*1.e-3
Xdist = Xdist-Xdist[0]
Xdist = Xdist/Xdist[-1]*Lsection*1.e-3  # normalized, km



import mod_colormaps as mclrmp
clrmp_name = 'winter'
clr_ramp   = [1, 1, 1]   # add white color at the end of the colormap
clrmp = mclrmp.addendclr_colormap(clrmp_name, clr_ramp, nramp=0.1, ramp_start=False)
clrmp.set_bad(color=[0.5, 0.5, 0.5])
rmin = -7000.
rmax = 0.

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
im1 = m.pcolormesh(xR, yR, AA, cmap=clrmp, vmin=rmin, vmax=rmax)

m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

# Plot section with distance markers:
m.plot(xSCT,ySCT,'r-',linewidth=2)
nn = len(Xdist)
dx = 200.
xpnt = 0.
for ii in range(nn):
  msz=3
  if Xdist[ii] >= xpnt:
    if xpnt == 0.:
      m.plot(xSCT[ii],ySCT[ii], marker='o', markersize=msz, color=(0.,0.,1))
    else:
      m.plot(xSCT[ii],ySCT[ii], marker='o', markersize=msz, color=(0.,0.,0.))
    xpnt = xpnt + dx

sttl = f"{sctnm} X markers dx={dx:4.0f} Max(X)={np.max(Xdist):5.0f}km"
ax1.set_title(sttl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
# extend: min, max, both
clb = plt.colorbar(im1, cax=ax2, orientation='vertical', extend='min')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.0f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

btx = 'plot_sections_ortho.py'
bottom_text(btx, fsz=6)


