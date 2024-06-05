"""
  Plot GOFS output T, S, or ssh 
  orthonormal projection
"""
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
#import mod_valid_utils as mvutil
importlib.reload(mcmp)

nrun   = "GOFS3.0"
nexpt  = 19.0
dnmb0  = mtime.datenum([1993,1,1])
varplt = 'srfhgt' # salin temp srfhgt 
kplt   = 0


with open('pypaths_gfdlpub.yaml') as ff:
  dct = yaml.safe_load(ff)

pthgofs = '/work/Dmitry.Dukhovskoy/GOFS3.0/expt_19.0/'

f_ssh   = False
f_temp  = False
f_salt  = False
if varplt == "srfhgt":
  f_ssh = True
elif varplt == "temp":
  f_temp = True
elif varplt == "salin":
  f_salt = True

# HYCOM GOFS grid:
pthgrid  = dct["GOFS3.0"]["19.0"]["pthgrid"]
ftopo    = dct["GOFS3.0"]["19.0"]["ftopo"]
fgrid    = dct["GOFS3.0"]["19.0"]["fgrid"]
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)
idmh = HH.shape[1]
jdmh = HH.shape[0]


# Read MOM6 NEP domain nx=342, ny=816
pthgrid_mom = dct["MOM6_NEP"]["test"]["pthgrid"]
ftopo_mom   = dct["MOM6_NEP"]["test"]["ftopo"]
fgrid_mom   = dct["MOM6_NEP"]["test"]["fgrid"]
pthoutp     = dct["MOM6_NEP"]["test"]["pthoutp"]
dftopo_mom  = os.path.join(pthgrid_mom, ftopo_mom)
dfgrid_mom  = os.path.join(pthgrid_mom, fgrid_mom)
LONM, LATM  = mom6util.read_mom6grid(dfgrid_mom)
HHM         = mom6util.read_mom6depth(dftopo_mom)
jdm         = np.shape(HHM)[0]
idm         = np.shape(HHM)[1]

DV      = mtime.datevec(dnmb0)
YR      = DV[0]
MM      = DV[1]
DD      = DV[2]
_, jday = mtime.dnmb2jday(dnmb0)
jday    = int(jday)
hr      = 15

# 190_archv.1993_001_15.a
fhycom = f"{nexpt*10:3.0f}_archv.{YR}_{jday:03d}_{hr:02d}"
fina   = os.path.join(pthgofs,fhycom) + '.a'
finb   = os.path.join(pthgofs,fhycom) + '.b'
A2d, _, _, _ = mhycom.read_hycom(fina, finb, varplt, rLayer=kplt+1)
A2d    = np.where(A2d > 1.e19, np.nan, A2d)
if f_ssh: A2d = A2d/9.806

# MOM domain bndry on ortho projection, need LON/LAT:
IBND = np.array([])
JBND = np.array([])
dstp = 10
# Western bndry
IBND = LONM[0:jdm:dstp, 0]
JBND = LATM[0:jdm:dstp, 0]
# northern bndry
IBND = np.append(IBND, LONM[jdm-1, 0:idm:dstp])
JBND = np.append(JBND, LATM[jdm-1, 0:idm:dstp])
# eastern bndry
IBND = np.append(IBND, LONM[0:jdm:dstp, idm-1])
JBND = np.append(JBND, LATM[0:jdm:dstp, idm-1])
# southern
IBND = np.append(IBND, LONM[0, 0:idm:dstp])
JBND = np.append(JBND, LATM[0, 0:idm:dstp])

# MOM domain bndry indices:
pthtmp  = '/work/Dmitry.Dukhovskoy/data/mom6_nep_tmp/'
fboutp  = 'mom_bndry_gofs.pkl'
dfboutp = os.path.join(pthtmp, fboutp)
print(f"Loading {dfboutp}")
with open(dfboutp, 'rb') as fid:
  IX, JX = pickle.load(fid)
imin = np.argmin(IX)
imax = np.argmax(IX)
jmin = np.argmin(JX)
jmax = np.argmax(JX)
Xv = np.array([IX[imin],IX[jmax],IX[imax],IX[jmin]])
Yv = np.array([JX[imin],JX[jmax],JX[imax],JX[jmin]])

if f_ssh:
  import mod_misc1 as mmisc
# Find regional mean and demean:
  IIH, JJH = np.meshgrid(np.arange(idmh), np.arange(jdmh))
  _, IP, JP = mmisc.inpolygon_v2(IIH, JJH, Xv, Yv)
  mnA = np.nanmean(A2d[JP,IP])
  A2d   = A2d - mnA
  cmpr  = mutil.colormap_ssh(nclrs=200)
  rmin = -0.5
  rmax = 0.5
  cntr1 = [x/10 for x in range(-10, 10, 1)]
elif f_salt:
  cmpr  = mcmp.colormap_haline()
  rmin  = 31.5
  rmax  = 35.5
  if kplt>15:
    rmin = 34.
    rmax = 34.9
  cntr1 = [x/10 for x in range(348, 358, 1)]
elif f_temp:
  cmpr = mcmp.colormap_temp_coldhot()
  rmin  = -2.
  rmax  = 25.
  if kplt>15:
    rmin = -1.
    rmax = 7.
  cntr1 = [x/10 for x in range(-10, 250, 10)]

nx     = A2d.shape[1]
ny     = A2d.shape[0]

ctitle = f'{nrun} {nexpt:3.1f} {varplt} Lr={kplt+1} {YR}/{MM:02d}/{DD:02d}'
cmpr.set_bad(color=[0.1, 0.1, 0.1])


# Set up orthographic projection
from mpl_toolkits.basemap import Basemap, cm

# Add extra row/col for plotting with pcolormesh
lonw = LON.copy()
latw = LAT.copy()
lonw = np.insert(lonw, -1, lonw[:,-1]+0.01, axis=1)
lonw = np.insert(lonw, -1, lonw[-1,:]+0.01, axis=0)
latw = np.insert(latw, -1, latw[-1,:]+0.01, axis=0)
latw = np.insert(latw, -1, latw[:,-1]+0.01, axis=1)

lon0 = 220.
lat0 = 50.
res  = 'l'
m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution=res)
xR, yR = m(lonw, latw)
PMsk = ( (xR > 1e20) | (yR > 1e20) )
AA = A2d.copy()
AA = np.insert(AA, -1, AA[:,-1], axis=1)
AA = np.insert(AA, -1, AA[-1,:], axis=0)
AA[PMsk] = np.nan
xR[PMsk]   = 1.e30
yR[PMsk]   = 1.e30

AA = AA[0:ny, 0:nx]

# MOM NEP domain:
xBND = []
yBND = []
xBND, yBND = m(IBND, JBND)

#================
#    Plotting

plt.ion()

fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
plt.clf()

ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
im1 = m.pcolormesh(xR, yR, AA, shading='flat', cmap=cmpr,\
                   vmin=rmin, vmax=rmax)

m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))
if len(xBND) > 0:
  m.plot(xBND, yBND, 'w.')


ax1.set_title(ctitle)


ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
             ax1.get_position().y0,0.02,
             ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, extend='both')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=10)

#fig1.colorbar(im,ax=ax1,orientation='horizontal')
btx = 'plotGOFS_momTS_ortho.py'
bottom_text(btx, pos=[0.02, 0.02])





