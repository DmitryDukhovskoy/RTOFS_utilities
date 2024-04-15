# Check horizontal interpolation 
# of GOFS onto  MOM 
#  Velocity fields
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
kk         = 1        # layer to plot 1, ..., kdmh
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

if fldint == "u-vel.":
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
#HH           = np.where(HH>=0, np.nan, HH)
jdh          = np.shape(HH)[0]
idh          = np.shape(HH)[1]

huge   = 1.e25
rg     = 9806.

# 190_archv.1993_001_15.a
fhycom = f"{nexpt*10:3.0f}_archv.{YR}_{jday:03d}_{hr:02d}"
fina   = os.path.join(pthrun,fhycom) + '.a'
finb   = os.path.join(pthrun,fhycom) + '.b'
U2d, idmh, jdmh, kdmh = mhycom.read_hycom(fina, finb, 'u-vel.', rLayer=kk)
U2d[np.where(U2d>huge)] = np.nan
V2d, idmh, jdmh, kdmh = mhycom.read_hycom(fina, finb, 'v-vel.', rLayer=kk)
V2d[np.where(V2d>huge)] = np.nan

# Read barotropic velocities:
# Only for archv output !
F, _, _, _  = mhycom.read_hycom(fina,finb,'u_btrop')
Ubtrop      = np.where(F>=huge, 0., F)
F, _, _, _  = mhycom.read_hycom(fina,finb,'v_btrop')
Vbtrop      = np.where(F>=huge, 0., F)

# Add depth-average U:
print(f"Adding depth-average U")
U2d = U2d + Ubtrop
V2d = V2d + Vbtrop

# Collocate U & V:
U2d = mhycom.collocateU2P(U2d)
V2d = mhycom.collocateV2P(V2d)
S2d = np.sqrt(U2d**2 + V2d**2)
U2d = np.where(HH>0., np.nan, U2d)
V2d = np.where(HH>0., np.nan, V2d)
S2d = np.where(HH>0., np.nan, S2d)

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
fluintrp  = f"gofs2mom_nep_hrzi-uvel_{grid_shape}_{rdate}.pkl"
dfluintrp = os.path.join(pthtmp,fluintrp)
flvintrp  = f"gofs2mom_nep_hrzi-vvel_{grid_shape}_{rdate}.pkl"
dflvintrp = os.path.join(pthtmp,flvintrp)

print(f"Loading {dfluintrp}")
with open(dfluintrp,'rb') as fid:
  U3di = pickle.load(fid)

print(f"Loading {dflvintrp}")
with open(dflvintrp,'rb') as fid:
  V3di = pickle.load(fid)

ikk  = kk -1
U2di = U3di[ikk,:,:].squeeze()
V2di = V3di[ikk,:,:].squeeze()
# Collocate MOM u & v:
U2di = mom6util.collocateU2H(U2di, grid_shape) 
V2di = mom6util.collocateV2H(V2di, grid_shape) 
S2di = np.sqrt(U2di**2 + V2di**2)
U2di = np.where(HHM>0., np.nan, U2di)
V2di = np.where(HHM>0., np.nan, V2di)
S2di = np.where(HHM>0., np.nan, S2di)

# Some checking
ih = 1405
jh = 2629
xh = LON[jh,ih]
yh = LAT[jh,ih]
im, jm = mutil.find_indx_lonlat(xh, yh, hlon, hlat)
print(f"HYCOM: u={U2d[jh,ih]} v={V2d[jh,ih]} |U|={S2d[jh,ih]}")
print(f"MOM6:  u={U2di[jm,im]} v={V2di[jm,im]} |U|={S2di[jm,im]}")

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
  if ikk>15:
    rmin = 34.
    rmax = 34.9
  cntr1 = [x/10 for x in range(348, 358, 1)]
elif fldint == 'temp':
#  cmpr = mcmp.colormap_temp()
  cmpr = mcmp.colormap_temp_coldhot()
  rmin  = -1.
  rmax  = 25.
  if ikk>15:
    rmin = -1.
    rmax = 7.
  cntr1 = [x/10 for x in range(-10, 250, 10)]
elif fldint == 'u-vel.' or fldint == 'v-vel.':
#  cmpr = mcmp.colormap_uv(nclrs=200)
#  rmin = -0.5
#  rmax = 0.5
  cmpr = mcmp.colormap_speed()
  rmin = 0.
  rmax = 0.5
  if ikk>15:
    rmin = 0
    rmax = 0.1

cmpr.set_bad(color=[0.9, 0.9, 0.9])

#================

plt.ion()

ctlt = f"Interp {fldint} GOFS --> MOM layer={kk} {rdate}"
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
im1 = ax1.pcolormesh(S2di, \
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

btx = 'check_UVinterp.py'
'check_hinterp.py'
bottom_text(btx)


# Plot hycom:
# Domain of interest:
xl1 = 1010
xl2 = 2280
yl1 = 1500
yl2 = 3050
ctlt2 = f"{nrun}-{expt} {fhycom} {fldint} layer {kk}"

fig2 = plt.figure(2, figsize=(9,8))
plt.clf()
ax21 = plt.axes([0.1, 0.1, 0.8, 0.8])
im2  = ax21.pcolormesh(S2d, \
                 cmap=cmpr,\
                 vmin=rmin, \
                 vmax=rmax)

ax21.axis('scaled')
ax21.set_xlim([xl1, xl2])
ax21.set_ylim([yl1, yl2])

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
clb = plt.colorbar(im2, cax=ax22, orientation='vertical', extend='both')
ax22.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax22.set_yticklabels(ax22.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

bottom_text(btx)


# ==================
f_vect = False
if f_vect:
  print('Plotting vectors')
  di  = 30
  Iv  = [ix for ix in range(0,idm,di)]
  Jv  = [ix for ix in range(0,jdm,di)]
  IJh = []
  IJm = []
  for iv in range(len(Iv)):
    for jv in range(len(Jv)):
      iv0 = Iv[iv]
      jv0 = Jv[jv]
      xmm = hlon[jv0,iv0]
      ymm = hlat[jv0,iv0]
      um  = U2di[jv0,iv0]
      vm  = V2di[jv0,iv0]
# Normalize:
      sm = np.sqrt(um**2 + vm**2)
      if sm <= 1.e-5: continue
      um = um/sm
      vm = vm/sm

      ax1.quiver(iv0, jv0, um, vm, units='inches', 
                 scale=2.5, width=0.025, color=[0.1,0.1,0.1])

      ihh, jhh = mutil.find_indx_lonlat(xmm, ymm, LON, LAT)
      uh = U2d[jhh,ihh]
      vh = V2d[jhh,ihh]
      sh = np.sqrt(uh**2 + vh**2) 
      uh = uh/sh
      vh = vh/sh
      ax21.quiver(ihh, jhh, uh, vh, units='inches',
                  scale=2.8, width=0.025, color=[0.1,0.1,0.1])
     

#      IJh.append([ihh,jhh])
#      IJm.append([iv0,jv0])








