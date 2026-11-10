"""
  Create ARC0.08 topo from GLBb0.08_09m11ob2_mom6.nc depth file
  !!! Note: Depth are positive in MOM6 (land = 0) !!!

    coordinates on "supergrid" (half grid points of actual grid)
    default - coordinates of h-points (scalar) are derived
    options: 'qgrid' - coordinated in the upper-right corner
                      q-points, vorticity points

   https://mom6.readthedocs.io/en/main/api/generated/pages/Discrete_Grids.html
 
   q(i-1,j)   ------------  v(i,j) --------------  q(i,j)
     |                                              |  
     |                                              |  
     |                                              |  
     |                                              |  
     |                                              |  

   u(i-1,j)               * h(i,j)                u(i,j)

     |                                              |  
     |                                              |  
     |                                              |  
     |                                              |  
   q(i-1,j-1) ------------  v(i,j-1) ------------  q(i,j-1)


  grid nonsymmetr: All arrays are declared with the same shape (isd:ied,jsd:jed)
                 There are fewer staggered variables to the south-west of the 
                 computational domain. An operator applied at h-point locations 
                 involving u- or v- point data can not have as wide a stencil on 
                 the south-west side of the processor domain as it can on 
                 the north-east side.

  grid symmetr:   Arrays have different shapes depending on their staggering 
                  location on the Arakawa C grid. 

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
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_misc1 as mmisc

pthin  = '/scratch4/NCEPDEV/stmp/Dmitry.Dukhovskoy/MOM6_CICE6_0.08global/fix'
pthout = '/scratch4/NCEPDEV/stmp/Dmitry.Dukhovskoy/ARC0.08/topo_grid'
f_save = True

# Read Global topo:
fltopo_in = 'depth_GLBb0.08_09m11ob2_mom6.nc'
dflin = os.path.join(pthin, fltopo_in)

print(f"Reading MOM6 topo: {dflin}")
with xarray.open_dataset(dflin) as ds:
  HH = ds["depth"].data

# h-grid lon/lats:
dfgrid = os.path.join(pthin, 'regional.mom6.nc')
LONh, LATh = mmom6.read_mom6grid(dfgrid, grdpnt='hgrid')

# Super-grid:
# Alfa: the counterclockwise angle from TRUE EAST to the model grid +x direction
with xarray.open_dataset(dfgrid) as dsgrid:
  LONS = dsgrid["x"].data.squeeze()  # degrees_east
  LATS = dsgrid["y"].data.squeeze()  # degrees north
  alfa = dsgrid["angle_dx"].data.squeeze()  # angle between model X axis and true East
  alfa_units = dsgrid["angle_dx"].attrs["units"] 
  dx   = dsgrid["dx"].data.squeeze()
  dy   = dsgrid["dy"].data.squeeze()
  area = dsgrid["area"].data.squeeze()
  tile = dsgrid["tile"].data 

lstring = 255

# Global dim
jdmG, idmG = HH.shape  
jdmS, idmS = LONS.shape

def subset_arc(HH, jcut, iL1, iL2, iR1, iR2):
  HL = HH[jcut:, iL1:iL2+1]
  HR = HH[jcut:, iR1:iR2+1]
  #print(f"HL shape: {HL.shape}")
  #print(f"HR shape: {HR.shape}")
  # Fold left half:
  HLf = np.flipud(np.fliplr(HL))
  HARC = np.concatenate((HR, HLf), axis=0)

  return HARC

def subset_arc_supergrid(HH, jcut, iL1, iL2, iR1, iR2, 
                         alfunits="degrees", angle=False, delete_seam=True):
  """
    Need to remove "seam" row in one of the 
    halves to avoid duplicating
    but only for arrays with values at all points on super grid
    along Y axis (e.g., dy - missing last values near the upper bndry)

    For true-east - North angle: rotate the upper half accordingly
           to match new grid X axis orientation
  """
  HL = HH[jcut:, iL1:iL2+1]
  HR = HH[jcut:, iR1:iR2+1]

  HLf = np.flipud(np.fliplr(HL))
  if angle:
    if alfunits == 'degrees':
      HLf = (HLf + 180.0) % 360.0
    else:
      #HLf = (HLf + np.pi) % (2*np.pi)
      #[-pi, pi]:
      HLf = (HLf + np.pi + np.pi) % (2*np.pi) - np.pi

  # remove duplicated seam row 
  if delete_seam:
    HARC = np.concatenate((HR, HLf[1:, :]), axis=0)
  else:
    HARC = np.concatenate((HR, HLf), axis=0)

  return HARC

# Set of points along the seam line in the right half of the grid:
# Left half flips up-down and goes "on top" of the right half in polar projection
imidR = idmG // 2
imidS = (idmS - 1) // 2


# Select domain:
# Indices in a regular grid h-points:
jcutBr = 2241
icutLr = 233  # cut from the mid-line
icutRr = 450  # cut from the L/R boundaries
iL1r = icutRr
iL2r = imidR - icutLr
iR1r = imidR + icutLr
iR2r = idmG - icutRr

LONh_sub = subset_arc(LONh, jcutBr, iL1r, iL2r, iR1r, iR2r)
LATh_sub = subset_arc(LATh, jcutBr, iL1r, iL2r, iR1r, iR2r)

# Convert to supergrid:
jcutB = jcutBr * 2
jcutU = LONS.shape[0]
iL1   = iL1r *2
iL2   = iL2r *2 + 2
iR1   = iR1r * 2 
iR2   = iR2r * 2 + 2

# Subset supergrid:
hlon = subset_arc_supergrid(LONS, jcutB, iL1, iL2, iR1, iR2)
hlat = subset_arc_supergrid(LATS, jcutB, iL1, iL2, iR1, iR2)

# Check
hlon_h = hlon[1::2, 1::2]
hlat_h = hlat[1::2, 1::2]
assert hlat_h.shape == LATh_sub.shape, "hlat shape error"
dlt = np.abs((hlon_h - LONh_sub + 180) % 360 - 180)
assert np.max(dlt) < 1.e-10, f"Subpergrid hlon subset h-points shifted"

# Check north pole seam at q-points:
jNP = hlon.shape[0] // 2 
DLT = np.sqrt(
    (hlon[jNP, :] - hlon[jNP-1, :])**2 +
    (hlat[jNP, :] - hlat[jNP-1, :])**2
)

print(f"Check: max hlon/hlat dlt at the seam expext <= 0.08: {np.nanmax(DLT)}")

# Compute dy at the NP seam that is missing:
#dy_seam = mmisc.dist_sphcrd(hlat[jNP, :],hlon[jNP, :],hlat[jNP-1, :],hlon[jNP-1, :])

# Subet other Grid vars:
dx_sub = subset_arc_supergrid(dx, jcutB, iL1, iL2-1, iR1, iR2-1) 
dy_sub = subset_arc_supergrid(dy, jcutB, iL1, iL2, iR1, iR2, delete_seam=False)
area_sub = subset_arc_supergrid(area, jcutB, iL1, iL2-1, iR1, iR2-1, delete_seam=False) 
alfa_sub = subset_arc_supergrid(alfa, jcutB, iL1, iL2, iR1, iR2, angle=True, alfunits=alfa_units)

# keep angle in [-180, 180]:
alfa_sub = (alfa_sub + 180.) % 360 - 180.

nyp, nxp = hlon.shape
ny = nyp - 1
nx = nxp - 1

assert hlat.shape == (nyp,nxp), f"hlat shape is wrong {hlat.shape}"
assert dx_sub.shape == (nyp,nx), f"dx_sub shape is wrong {dx_sub.shape}"
assert dy_sub.shape == (ny,nxp), f"dy_sub shape is wrong {dy_sub.shape}"
assert area_sub.shape == (ny,nx), f"area_sub shape is wrong {area_sub.shape}"


if f_save:
  fltopo_out = 'regional.ARC008_09m11ob2.mom6.nc'
  dflout = os.path.join(pthout, fltopo_out)

  # need to convert "tile1" into a char array of length lstring
  tile_str = tile.tobytes().decode('utf-8').strip()
  tile_arr = np.array(list(tile_str.ljust(lstring)), dtype='S1')

  ds_out = xarray.Dataset(
    data_vars=dict(
      x=(["nyp", "nxp"], hlon.astype(np.float32)),
      y=(["nyp", "nxp"], hlat.astype(np.float32)),
      angle_dx=(["nyp","nxp"], alfa_sub.astype(np.float32)),
      dx=(["nyp","nx"], dx_sub.astype(np.float32)),
      dy=(["ny","nxp"], dy_sub.astype(np.float32)),
      area=(["ny","nx"], area_sub.astype(np.float32)),
      tile=(["string"], tile_arr), 
      )
  )

  # Attributes
  ds_out["x"].attrs = {
      "units": "degrees_east",
  }

  ds_out["y"].attrs = {
      "units": "degrees_north",
  }

  ds_out["angle_dx"].attrs = {
      "units": "degrees",
  }

  ds_out["dy"].attrs = {
      "units": "m",
  }

  ds_out["dx"].attrs = {
      "units": "m",
  }

  ds_out["area"].attrs = {
      "units": "m2",
  }

  # Global attributes:
  from datetime import datetime
  ds_out.attrs.update({
      "title": f"ARC0.08 regional grid from global {fltopo_in}",
      "source": "create_coordARC008.py",
      "history": f"Modified {datetime.now().isoformat()}",
  })

  encoding = {
      "x": {"dtype": "float32"},
      "y": {"dtype": "float32"},
      "angle_dx": {"dtype": "float32"},
      "dx": {"dtype": "float32"},
      "dy": {"dtype": "float32"},
      "area": {"dtype": "float32"},
  }

  print(f"Saving ARC topo --> {dflout}")
  ds_out.to_netcdf(dflout, encoding=encoding)
  


check_plt = False
if check_plt:
  plt.ion()

  fplt = 'alfa'
  fplt = 'area'
  fplt = 'dx'
  fplt = 'dy'
  if fplt == 'alfa':
    AA = alfa_sub.copy()
    clrmp = mclrmps.colormap_uv()
    rmin = -180.
    rmax = 180.
  elif fplt == 'area':
    AA = area_sub.copy() * 1e-6 # m2 --> km2
    clrmp = mclrmps.colormap_haline2()
    rmin = 0.
    rmax = 9.
  elif fplt == 'dx':
    AA = dx_sub.copy()  # m
    clrmp = mclrmps.colormap_haline2()
    rmin = 1000.
    rmax = 3000.
  elif fplt == 'dy':
    AA = dy_sub.copy()  # m
    clrmp = mclrmps.colormap_haline2()
    rmin = 1000.
    rmax = 3000.

  jdm, idm = AA.shape
 
  clrmp.set_bad(color=[0.1, 0.1, 0.1])

  fig1 = plt.figure(1,figsize=(9,9))
  plt.clf()
  ax1 = plt.axes([0.08, 0.08, 0.8, 0.8])

  ax1.cla()
  img = ax1.pcolormesh(AA, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')
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
      levels=range(-360, 359, 10),
      colors=[(0.8, 0.8, 0.8)],
      linewidths=0.5,
      linestyles='solid'
  )

  sttl = f"ARC0.08 Super-grid, {fplt} {jdm} x {idm},\nmin={np.min(AA):.3f}, max={np.max(AA):.3f}"
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


