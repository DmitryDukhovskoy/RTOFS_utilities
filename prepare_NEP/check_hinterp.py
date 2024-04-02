# Check horizontal interpolation 
# of GOFS onto  MOM 
# 
# Plot NEP domain on GOFS grid
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
#import struct
import datetime
import pickle
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import time
import timeit
import yaml
from netCDF4 import Dataset as ncFile

#PPTHN = '/home/Dmitry.Dukhovskoy/python'
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

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_read_hycom as mhycom
import mod_colormaps as mcmp
import mod_mom6 as mom6util
import mod_regmom as mrgm
#import mod_valid_utils as mvutil
importlib.reload(mcmp)

nrun       = "GOFS3.0"
nexpt      = 19.0
expt       = f"{nexpt:2.1f}"
YR         = 1993
jday       = 1
hr         = 15
fldint     = "u-vel."  # temp salin u-vel. v-vel. thknss
grid_shape = 'symmetr'   # MOM grid: symmetr/nonsymmetr

dnmb    = mtime.jday2dnmb(YR, jday)
dv      = mtime.datevec(dnmb)
ds      = mtime.datestr(dnmb)
MM      = dv[1]
DM      = dv[2]
rdate   = f"{YR}{MM:02d}{DM:02d}"

print(f"Plotting {fldint}")

pthout  = '/work/Dmitry.Dukhovskoy/data/mom6_nep_restart/'
pthtmp  = '/work/Dmitry.Dukhovskoy/data/mom6_nep_tmp/'
pthgofs = '/work/Dmitry.Dukhovskoy/data/GOFS3.0/expt_19.0/'

if fldint == "temp":
  grid_var = 'hgrid'
  fldiout  = 'temp'
elif fldint == "salin":
  grid_var = 'hgrid'
  fldiout  = 'salin'
elif fldint == "u-vel.":
  grid_var = 'ugrid'
  fldiout  = 'uvel'
elif fldint == "v-vel.":
  grid_var = 'vgrid'
  fldiout  = 'vvel'
elif fldint == "thknss":
  grid_var = 'hgrid'
  fldiout  = 'lrthknss'


with open('pypaths_gfdlpub.yaml') as ff:
  dct = yaml.safe_load(ff)

pthrun  = dct[nrun][expt]["pthrun"]
pthgrid = dct[nrun][expt]["pthgrid"]
ftopo   = dct[nrun][expt]["ftopo"]
fgrid   = dct[nrun][expt]["fgrid"]

LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
HH           = np.where(HH>=0, np.nan, HH)
jdh          = np.shape(HH)[0]
idh          = np.shape(HH)[1]

huge   = 1.e25
rg     = 9806.

# 190_archv.1993_001_15.a
fhycom = f"{nexpt*10:3.0f}_archv.{YR}_{jday:03d}_{hr:02d}"
fina   = os.path.join(pthrun,fhycom) + '.a'
finb   = os.path.join(pthrun,fhycom) + '.b'
A3d, idmh, jdmh, kdmh = mhycom.read_hycom(fina,finb,fldint)
A3d[np.where(A3d>huge)] = np.nan

#
# Read MOM6 NEP domain nx=342, ny=816
pthgrid_mom = dct["MOM6_NEP"]["test"]["pthgrid"]
ftopo_mom   = dct["MOM6_NEP"]["test"]["ftopo"]
fgrid_mom   = dct["MOM6_NEP"]["test"]["fgrid"]
dftopo_mom  = os.path.join(pthgrid_mom, ftopo_mom) 
dfgrid_mom  = os.path.join(pthgrid_mom, fgrid_mom) 
hlon, hlat  = mom6util.read_mom6grid(dfgrid_mom, grid=grid_shape, grdpnt=grid_var)
HHM         = mom6util.read_mom6depth(dftopo_mom) 
jdm         = np.shape(HHM)[0]
idm         = np.shape(HHM)[1]
# Convert to -180/180:
hlon   = np.where(hlon > 180.0, hlon-360., hlon)

# Output interpolated fields on hycom layers:
flintrp  = f"gofs2mom_nep_hrzi-{fldiout}_{grid_shape}_{rdate}.pkl"
dflintrp = os.path.join(pthtmp,flintrp)

print(f"Loading {dflintrp}")
with open(dflintrp,'rb') as fid:
  A3di = pickle.load(fid)


kk = 10
AA = A3di[kk,:,:].squeeze()
AH = A3d[kk,:,:].squeeze()

# MOM domain bndry:
fboutp  = 'mom_bndry_gofs.pkl'
dfboutp = os.path.join(pthtmp, fboutp)
print(f"Loading {dfboutp}")
with open(dfboutp, 'rb') as fid:
  IBND, JBND = pickle.load(fid)


import mod_colormaps as mcmp
importlib.reload(mcmp)

if fldint == 'salin':
#  cmpr = mcmp.colormap_salin(clr_ramp=[0.94,0.95,1])
  cmpr  = mcmp.colormap_haline()
  rmin  = 31.5
  rmax  = 35.5
  if kk>15:
    rmin = 34.
    rmax = 34.9
  cntr1 = [x/10 for x in range(348, 358, 1)]
elif fldint == 'temp':
#  cmpr = mcmp.colormap_temp()
  cmpr = mcmp.colormap_temp_coldhot()
  rmin  = -1.
  rmax  = 25.
  if kk>15:
    rmin = -1.
    rmax = 7.
  cntr1 = [x/10 for x in range(-10, 250, 10)]
elif fldint == 'u-vel.' or fldint == 'v-vel.':
  cmpr = mutil.colormap_ssh(nclrs=200)
  rmin = -0.5
  rmax = 0.5
  if kk>15:
    rmin = -0.2
    rmax = 0.2

cmpr.set_bad(color=[0.1, 0.1, 0.1])

#================

plt.ion()

ctlt = f"Interp {fldint} GOFS --> MOM layer={kk+1} {rdate}"
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
im1 = ax1.pcolormesh(AA, \
                 cmap=cmpr,\
                 vmin=rmin, \
                 vmax=rmax)

ax1.axis('scaled')
ax1.set_xlim([0, idm-1])
ax1.set_ylim([0, jdm-1])

LON1 = hlon.copy()
LON1 = np.where(LON1<0., hlon+360., hlon)
clrg = [(0.95,0.95,0.95)]
ax1.contour(LON1,list(range(100,320,10)),
          colors=clrg,
          linestyles='solid',
          linewidths=1.0)
ax1.contour(hlat,list(range(0,89,10)),
          colors=clrg,
          linestyles='solid',
          linewidths=1.0)

ax1.set_title(ctlt)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, orientation='vertical', extend='both')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

btx = 'check_hinterp.py'
bottom_text(btx)


# Plot hycom:
ctlt2 = f"{nrun}-{expt} {fhycom} {fldint} layer {kk+1}"

fig2 = plt.figure(2, figsize=(9,8))
plt.clf()
ax21 = plt.axes([0.1, 0.1, 0.8, 0.8])
im2  = ax21.pcolormesh(AH, \
                 cmap=cmpr,\
                 vmin=rmin, \
                 vmax=rmax)

ax21.axis('scaled')
ax21.set_xlim([1010, 2280])
ax21.set_ylim([1500, 3050])

ax21.plot(IBND, JBND, 'w.')
ax21.set_title(ctlt2)

LON1 = LON.copy()
LON1 = np.where(LON1<0., LON+360., LON)
clrg = [(0.95,0.95,0.95)]
ax21.contour(LON1,list(range(100,181,10)),
          colors=clrg,
          linestyles='solid',
          linewidths=1.0)
ax21.contour(LON,list(range(-170,-60,10)),
          colors=clrg,
          linestyles='solid',
          linewidths=1.0)
ax21.contour(LAT,list(range(0,89,10)),
          colors=clrg,
          linestyles='solid',
          linewidths=1.0)

ax22 = fig2.add_axes([ax21.get_position().x1+0.025, ax21.get_position().y0,
                   0.02, ax21.get_position().height])
clb = plt.colorbar(im1, cax=ax22, orientation='vertical', extend='both')
ax22.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax22.set_yticklabels(ax22.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

bottom_text(btx)








