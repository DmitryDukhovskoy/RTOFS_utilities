"""
  Find point along the N. Pole seam that corresponds to given pnt on the seam
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
from mpl_toolkits.basemap import Basemap, cm
import argparse

PPTHN = None
if 'PPTHN' not in locals() or PPTHN is None:
  cwd = os.getcwd()
  parts = cwd.split(os.sep)
  if 'python' in parts:
    idx = parts.index('python')
    PPTHN = os.sep + os.path.join(*parts[:idx + 1])
  else:
    raise RuntimeError("Directory 'python' not found in current working directory path.")

sys.path.extend([
    os.path.join(PPTHN, 'MyPython', 'hycom_utils'),
    os.path.join(PPTHN, 'MyPython', 'draw_map'),
    os.path.join(PPTHN, 'MyPython'),
    os.path.join(PPTHN, 'MyPython', 'mom6_utils')
])

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_read_hycom as mhycom

init_date = 20250704  #
init_hr = 0
regn = 'north'
fhr = 0


# Init date:
dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
yrI, mmI, ddI, hrI = mtime.datevec(dnmbI)[:4]


pthice = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/RTOFS/cice4_fcast/{init_date}"
if fhr == 0:
  flice = f"rtofs_glo.t{init_hr:02d}z.n00.cice_inst.nc"
else:
  flice = f"rtofs_glo.t{init_hr:02d}z.f{fhr:02d}.cice_inst.nc"
dflice = os.path.join(pthice, flice)

# Get grid
with xarray.open_dataset(dflice) as dcice:
  LON  = dcice["TLON"].data
  LAT  = dcice["TLAT"].data
  hice = dcice["hi"].data.squeeze()

JDIM, IDIM = LON.shape
JDIM = JDIM + 1   # ocean grid has + 1 row

# Read RTOFS topo:
# Note that RTOFS grid has +1 row at the top compared to CICE6
pthtopo = '/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/RTOFS/topo_grid/'
ftopo  = 'depth_GLBb0.08_09m11'
HH0 = mhycom.read_topo(pthtopo, ftopo, IDIM, JDIM)
HH0 = HH0[:-1,:]     # discard the extra row

# Subset hemispheres:
if regn == 'north':
  ilat = np.argmax(np.any(LAT >= 50, axis=1))
  hlat = LAT[ilat:,:]
  hlon = LON[ilat:,:]
  A2d  = hice[ilat:,:]
  HH   = HH0[ilat:,:]

elif regn == 'south':
  ilat = np.argmax(np.any(LAT >= -50, axis=1))
  hlat = LAT[:ilat,:]
  hlon = LON[:ilat,:]
  A2d  = hice[:ilat,:]
  HH   = HH0[ilat:,:]

jdm, idm = hlon.shape



i0 = 3312
j0 = jdm-1

# Set of points along the seam line in the right half of the grid:
imid = idm // 2
IR = np.arange(imid, imid + 2000, 20)
JR = np.zeros_like(IR) + j0

# Find indices in the left half of the grid corresponding to the RH indices:
IL = imid - (IR - imid + 1)
JL = np.zeros_like(IL) + j0


# Check segments across the discontinuity:
isctR = 3190
JSR = np.arange(1060, jdm)
ISR = np.zeros_like(JSR) + isctR

# Find corresponding indices in the Left part of the grid:
isctL = imid - (isctR - imid + 1)
JSL = np.flipud(JSR)
ISL = np.zeros_like(JSL) + isctL

clrmp = mclrmps.colormap_ice_thkn()
rmin = 0.
rmax = 3.
clrmp.set_bad(color=[0.1, 0.1, 0.1])
cntr_clr = [0.9,0.,1]



plt.ion()

if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)
elif regn == 'north':
  m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
  parallels = np.arange(50, 90, 5)
  meridians = np.arange(-360, 359., 45.)

xh, yh = m(hlon,hlat) # GFS coords


print("Plotting ...")

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.08, 0.1, 0.8, 0.8])

# Index space:
#ax1.axis('equal')
ax1.pcolormesh(A2d, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')
ax1.set_aspect('equal')
plt_sct = True
if plt_sct:
  ax1.plot(ISR, JSR, '.', color=[0.5,0.5,0.5])
  ax1.plot(ISL, JSL, '.', color=[0.,0.9,0.8])

plt_bndry = False
if plt_bndry:
  ax1.plot(IR, JR, '.', color=[0.5,0.5,0.5])
  ax1.plot(IL, JL, '.', color=[0.,0.9,0.8])

#  Spherical projection:
plt.clf()
ax1 = plt.axes([0.08, 0.1, 0.8, 0.8])


#m.drawcoastlines()
m.drawparallels(parallels,labels=[0,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10)

img = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')
#ax1.contour(xh, yh, HH, [0], linestyles='solid', colors=[cntr_clr], linewidths=1)

if plt_sct:
  ax1.plot(xh[JSR,ISR], yh[JSR,ISR], '.', color=[0.5,0.5,0.5])
  ax1.plot(xh[JSL,ISL], yh[JSL,ISL], '.', color=[0.,0.9,0.8])


if plt_bndry:
  ax1.plot(xh[JR,IR], yh[JR,IR], '.', color=[0.5,0.5,0.5])
  ax1.plot(xh[JL,IL], yh[JL,IL], '.', color=[0.,0.9,0.8])


ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
if rmin < 0:
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
else:
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')

ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=14)
clb.ax.tick_params(direction='in', length=12)

btx = 'find_RTOFS_polepnt.py'
bottom_text(btx, pos=[0.2, 0.01])



# Get x-section:
Asct = A2d[JSR,ISR]
Asct = np.append(Asct, A2d[JSL,ISL])
sttl = f'ithkn, RTOFS, init: {init_date}, fhr={fhr}, xsect iR={isctR} iL={isctL}'
ax1.set_title(sttl)


