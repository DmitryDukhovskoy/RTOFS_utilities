"""
  Correct "seam problem" along the grid N. Boundary
  in CICE4 fields
  For interpolation onto mesh025 grid
  Fields: 
  hi   - grid cell mean ice thickn
  hs   - -"-  -"- -"-  snow thickn
  aice - aggregated ice area
  Tsfc - snow/ice surf T 

  RTOFS CICE4 forecasts

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

init_date = 20250704  #
init_hr = 0
fhr = 0

parser = argparse.ArgumentParser()
parser.add_argument("--init", help=f"init date", choices=[20251231, 20250704], required=True, type=int)
parser.add_argument("--ihr", help=f"init hour, default {init_hr}", type=int)
args = parser.parse_args()

init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr


# Init date:
dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
yrI, mmI, ddI, hrI = mtime.datevec(dnmbI)[:4]


pthice = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/RTOFS/cice4_fcast/{init_date}"
if fhr == 0:
  flice = f"rtofs_glo.t{init_hr:02d}z.n00.cice_inst.nc"
else:
  flice = f"rtofs_glo.t{init_hr:02d}z.f{fhr:02d}.cice_inst.nc"
dflice = os.path.join(pthice, flice)

print(f"Reading restart: {dflice}")
ds_in = xarray.open_dataset(dflice)
ds_out = ds_in.copy(deep=True)
ds_in.close()

LON   = ds_out["TLON"].data
LAT   = ds_out["TLAT"].data
hice  = ds_out["hi"].data.squeeze()
hsnow = ds_out["hs"].data.squeeze()
aice  = ds_out["aice"].data.squeeze()

hice_new  = hice.copy()
hsnow_new = hsnow.copy()
aice_new  = aice.copy() 

jdm, idm = LON.shape

# Read RTOFS topo:
# Note that RTOFS grid has +1 row at the top compared to CICE6
pthtopo = '/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/RTOFS/topo_grid/'
ftopo  = 'depth_GLBb0.08_09m11'
HH = mhycom.read_topo(pthtopo, ftopo, idm, jdm+1)
HH = HH[:-1,:]     # discard the extra row
#LMsk = np.where(Himid + (H0>=0, 0, 1)

# Date :
dnmbP = dnmbI + fhr // 24
YR, MM, DD = mtime.datevec(dnmbP)[:3]

def interp_fld (A2d, JL, IL, JR, IR, idm, jdm, Xp, XI, wt_bary):
  Yp = np.concatenate([A2d[JL, IL], A2d[JR, IR]])

  Yint = []
  for xx in XI: 
    Pn = mintrp.barycentr_lagr_polynom(Xp, Yp, xx, wt=wt_bary) 
    Yint.append(Pn)

  return np.array(Yint)

# Set of points along the seam line in the right half of the grid:
imid = idm // 2

print("Start interpolation ...")
wt_bary = None
dlty = 6     # how many points from the N bndry seam
nfix = 3     # how many points to fix from the seam
jj0 = jdm - 1
for il in range(imid):
  if HH[-1,il] >= 0:
    continue

  JL = np.arange(jj0-dlty, jdm-nfix)
  IL = np.zeros_like(JL) + il

  ir = 2*imid - il - 1
  JR = np.flipud(JL)
  IR = np.zeros_like(JR) + ir

  # Distances along the sections:
  DL = JL - jdm
  DR = jdm - JR - 1

  Xp = np.concatenate([DL,DR])

  # Interpolation points:
  XI = np.arange(JL[-1]-jdm+1, jdm-JR[0]-1)  

  if wt_bary is None:
    wt_bary = mintrp.bary_weights(Xp)

  Yhi = interp_fld(hice,  JL, IL, JR, IR, idm, jdm, Xp, XI, wt_bary) 
  Yai = interp_fld(aice,  JL, IL, JR, IR, idm, jdm, Xp, XI, wt_bary) 
  Yhs = interp_fld(hsnow, JL, IL, JR, IR, idm, jdm, Xp, XI, wt_bary) 
  
  hice_new[-nfix:,il] = Yhi[:3]
  hice_new[-nfix:,ir] = np.flipud(Yhi)[:3]

  aice_new[-nfix:,il] = Yai[:3]
  aice_new[-nfix:,ir] = np.flipud(Yai)[:3]

  hsnow_new[-nfix:,il] = Yhs[:3]
  hsnow_new[-nfix:,ir] = np.flipud(Yhs)[:3]

hice_new = np.expand_dims(hice_new, axis=0)
hsnow_new = np.expand_dims(hsnow_new, axis=0)
aice_new = np.expand_dims(aice_new, axis=0)

ds_out["hi"].values[:]   = hice_new
ds_out["hs"].values[:]   = hsnow_new
ds_out["aice"].values[:] = aice_new

assert ds_out["hi"].shape == hice_new.shape, "Check shape of hi "
assert ds_out["hs"].shape == hsnow_new.shape, "Check shape of hs "
assert ds_out["aice"].shape == aice_new.shape, "Check shape of aice "

# Attributes:
from datetime import datetime
ds_out.attrs.update({
    "title": f"RTOFS CICE4 corrected ice & snow discontinuity across north boundary",
    "source": "correct_RTOFS_NPole_seam.py",
    "info": "RTOFS init_date {initdate} {flice}",
    "history": f"Modified {datetime.now().isoformat()}",
})

# Save:
base = os.path.splitext(flice)[0]
flrst_out = f"{base}.seam_crct.nc"

pth_out = pthice
dflrst_out = os.path.join(pth_out,flrst_out)
print(f"Saving RTOFS CICE4  --> {dflrst_out}")
ds_out.to_netcdf(dflrst_out, encoding={var: {'_FillValue': None} for var in ds_out.data_vars}, format='NETCDF3_64BIT')
ds_out.close()




check_plt = False
if check_plt:
  plt.ion()

  il = 1450
  jj0 = jdm - 1
  JLL = np.arange(jj0-dlty, jdm).astype(int)
  ILL = np.zeros_like(JLL) + il

  ir = 2*imid - il - 1
  JRR = np.flipud(JLL)
  IRR = np.zeros_like(JRR) + ir

  # Distances along the sections:
  DLL = JLL - jdm
  DRR = jdm - JRR - 1

  IX = np.concatenate([DLL,DRR])

  fig1 = plt.figure(1,figsize=(9,9))
  plt.clf()

  AA = hice
  AN = hice_new
  Yold = np.concatenate([AA[JLL,ILL], AA[JRR,IRR]])
  Ynew = np.concatenate([AN[JLL,ILL], AN[JRR,IRR]])

  ax1 = plt.axes([0.08, 0.7, 0.8, 0.25])
  sttl = f"hice, iL={il} - iR={ir}, baryc. Lagr. interp. polynom."
  ax1.plot(IX, Yold, '.-')
  ax1.plot(IX, Ynew, '.-')
  ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
  ax1.set_title(sttl) 

  AA = aice
  AN = aice_new
  Yold = np.concatenate([AA[JLL,ILL], AA[JRR,IRR]])
  Ynew = np.concatenate([AN[JLL,ILL], AN[JRR,IRR]])

  ax2 = plt.axes([0.08, 0.39, 0.8, 0.25])
  sttl = f"aice, iL={il} - iR={ir}, baryc. Lagr. interp. polynom."
  ax2.plot(IX, Yold, '.-')
  ax2.plot(IX, Ynew, '.-')
  ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
  ax2.set_title(sttl) 

  AA = hsnow
  AN = hsnow_new
  Yold = np.concatenate([AA[JLL,ILL], AA[JRR,IRR]])
  Ynew = np.concatenate([AN[JLL,ILL], AN[JRR,IRR]])

  ax3 = plt.axes([0.08, 0.08, 0.8, 0.25])
  sttl = f"hsnow, iL={il} - iR={ir}, baryc. Lagr. interp. polynom."
  ax3.plot(IX, Yold, '.-')
  ax3.plot(IX, Ynew, '.-')
  ax3.grid(True, which='both', linestyle='--', linewidth=0.5)
  ax3.set_title(sttl) 

  btx = 'correct_RTOFS_NPole_seam.py'
  bottom_text(btx, pos=[0.1, 0.04])


