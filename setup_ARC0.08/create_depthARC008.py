"""
  Create ARC0.08 topo from GLBb0.08_09m11ob2_mom6.nc depth file
  !!! Note: Depth are positive in MOM6 (land = 0) !!!

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
import mod_read_hycom as mhycom
import mod_interp1D as mintrp
import mod_mom6 as mmom6

pthin  = '/scratch4/NCEPDEV/stmp/Dmitry.Dukhovskoy/MOM6_CICE6_0.08global/fix'
pthout = '/scratch4/NCEPDEV/stmp/Dmitry.Dukhovskoy/ARC0.08/topo_grid'
f_save = True

# Read Global topo:
fltopo_in = 'depth_GLBb0.08_09m11ob2_mom6.nc'
dflin = os.path.join(pthin, fltopo_in)

print(f"Reading MOM6 topo: {dflin}")
with xarray.open_dataset(dflin) as ds:
  HH = ds["depth"].data

dfgrid = os.path.join(pthin, 'regional.mom6.nc')
LON, LAT = mmom6.read_mom6grid(dfgrid, grdpnt='hgrid')


# Global dim
jdmG, idmG = HH.shape  


def subset_arc(HH, jcut, iL1, iL2, iR1, iR2):
  HL = HH[jcut:, iL1:iL2+1]
  HR = HH[jcut:, iR1:iR2+1]
  #print(f"HL shape: {HL.shape}")
  #print(f"HR shape: {HR.shape}")
  HARC = np.concatenate((HR,np.flipud(np.fliplr(HL))), axis=0)

  return HARC

# Set of points along the seam line in the right half of the grid:
# Left half flips up-down and goes "on top" of the right half in polar projection
imid = idmG // 2

# Select domain:
jcut = 2241
icutL = 233  # cut from the mid-line
icutR = 450  # cut from the L/R boundaries
iL1 = icutR
iL2 = imid - icutL
iR1 = imid + icutL
iR2 = idmG - icutR

#HL = HH[jcut:, :imid]
#HR = HH[jcut:, imid:]
HARC = subset_arc(HH, jcut, iL1, iL2, iR1, iR2)
HPLT = -HARC.copy()
HPLT[HPLT>=0] = np.nan
jdm, idm = HPLT.shape

hmin = np.min(HARC[HARC>0])
hmax = np.max(HARC)
hland = np.min(HARC)
print(f"ARC0.08 shape: {jdm} x {idm}")
print(f"Max depth: {hmax:.4f}m, min depth: {hmin:.4f}m, land: {hland:.4f}m")

# Grid:
hlon = subset_arc(LON, jcut, iL1, iL2, iR1, iR2)
hlat = subset_arc(LAT, jcut, iL1, iL2, iR1, iR2)

wet_ARC = np.where(HARC<0, 1, 0)

if f_save:
  fltopo_out = 'depth_ARC0.08_09m11ob2_mom6.nc'
  dflout = os.path.join(pthout, fltopo_out)

  ds_out = xarray.Dataset(
      data_vars=dict(
          depth=(["ny", "nx"], HARC.astype(np.float32)),
          wet=(["ny", "nx"], wet_ARC.astype(np.float32)),
      )
  )

 # Define dimensions:
  ds_out = ds_out.assign_coords(
      ny=np.arange(HARC.shape[0]),
      nx=np.arange(HARC.shape[1])
  )

  # Attributes
  ds_out["depth"].attrs = {
      "units": "m",
      "description": "topographic depth at T-cell centers"
  }

  ds_out["wet"].attrs = {
      "units": "m",
      "description": "Values: 1=Ocean, 0=Land"
  }

  # Global attributes:
  from datetime import datetime
  ds_out.attrs.update({
      "title": f"ARC0.08 depth from global {fltopo_in}",
      "source": "create_depthARC008.py",
      "history": f"Modified {datetime.now().isoformat()}",
  })

  print(f"Saving ARC topo --> {dflout}")
  ds_out.to_netcdf(dflout)
  


check_plt = False
if check_plt:
  plt.ion()

  clrmp = mclrmps.colormap_topo1()
  rmin = -6000
  rmax = 0
  clrmp.set_bad(color=[0.1, 0.1, 0.1])

  fig1 = plt.figure(1,figsize=(9,9))
  plt.clf()
  ax1 = plt.axes([0.08, 0.08, 0.8, 0.8])
  ax1.contour(HH,[0], colors=[(0.7,0.7,0.7)])
  ax1.contour(LAT,[x for x in range(30,85,5)], colors=[(0., 0.5, 0.9)])
  ax1.plot([0, idmG],[jcut,jcut],'r-')
  ax1.plot([iL1, iL1],[jcut, jdmG],'r-')
  ax1.plot([iL2, iL2],[jcut, jdmG],'r-')
  ax1.plot([iR1, iR1],[jcut, jdmG],'r-')
  ax1.plot([iR2, iR2],[jcut, jdmG],'r-')
  ax1.plot([imid, imid],[0, jdmG],'k--')

  ax1.cla()
  img = ax1.pcolormesh(HPLT, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')
  ax1.axis('scaled')
  ax1.contour(
      hlat,
      levels=range(40, 88, 5),
      colors=[(0.8, 0.8, 0.8)],
      linewidths=0.5,
      linestyles='solid'
  )
  ax1.contour(
      hlon,
      levels=range(-170, 177, 10),
      colors=[(0.8, 0.8, 0.8)],
      linewidths=0.5,
      linestyles='solid'
  )

  sttl = f"ARC0.08 topo, dim: {jdm} x {idm}\nHmin={hmin:.3f}m, Hmax={hmax:.3f}, Land={hland:.3f}m"
  ax1.set_title(sttl)


  ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                     0.02, ax1.get_position().height])
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='min')
  
  ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=14)
  clb.ax.tick_params(direction='in', length=12)


  btx = 'create_depthARC008.py'
  bottom_text(btx, pos=[0.15, 0.02])


