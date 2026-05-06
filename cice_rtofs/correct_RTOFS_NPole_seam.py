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
import mod_mom6 as mmom

init_date = 20250704  #
init_hr = 0
fhr = 0

parser = argparse.ArgumentParser()
parser.add_argument("--init", help=f"init date", choices=[20251231, 20250704], required=True, type=int)
parser.add_argument("--ihr", help=f"init hour, default {init_hr}", type=int)
parser.add_argument("--fhr", help=f"f/cast hour, default=0", 
                   choices=[0,24,48,72,96,120,144,168,192], default=0,  type=int)
parser.add_argument("--save", help="Save corrected fields or not (for debugging or plotting), default=1",
           choices=[0,1], default=1, type=int)
args = parser.parse_args()

init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
fhr       = args.fhr
fsave = args.save == 1

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
LON   = ds_in["TLON"].data
LAT   = ds_in["TLAT"].data
hice  = ds_in["hi"].data.squeeze()
hsnow = ds_in["hs"].data.squeeze()
aice  = ds_in["aice"].data.squeeze()
ds_out = ds_in.copy(deep=True)
ds_in.close()

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

# Date :
dnmbP = dnmbI + fhr // 24
YR, MM, DD = mtime.datevec(dnmbP)[:3]

def interp_fld (A2d, JL, IL, JR, IR, idm, jdm, Xp, XI, wt_bary):
  Yp = np.concatenate([A2d[JL, IL], A2d[JR, IR]])
  ILR = np.concatenate([IL, IR])
  JLR = np.concatenate([JL, JR])

  dij = 20
  if np.any(np.isnan(Yp)):
    # substitute nans with means:
    inan = np.where(np.isnan(Yp))[0]
    for ii in inan:
      j0 = JLR[ii]
      i0 = ILR[ii]
      i1 = np.max([i0 - dij, 0])
      i2 = np.min([i0 + dij, idm])
      j1 = np.max([j0 - dij, 0])
      j2 = np.min([j0 + dij, jdm])
      Amean = np.nanmean(A2d[j1:j2, i1:i2])
      if np.isnan(Amean): Amean = 0.
      Yp[ii] = Amean
      
  Yint = []
  for xx in XI: 
    Pn = mintrp.barycentr_lagr_polynom(Xp, Yp, xx, wt=wt_bary) 
    Yint.append(Pn)

  if np.any(np.isnan(Yint)):
    print(f"WARN: Barycentr. Lagr. interpolation: NaNs in interpolated values")
  
  Yint = np.array(Yint)

  # No negative values
  Yint[Yint < 0.] = 0.

  return Yint

def solve_laplace(A2d, imid, dlty, nfix, Niter, eps_stop=0.01):
  """
    Iteratively solve grad^2(U) = 0
    stop when change over the gap region does not change > eps_stop

    Laplace smoothing over the gap region 
    - Ocean/valid data: fixed (Dirichlet BC)
    - Gap region: updated by Laplace iteration
    - Land (NaN mask): excluded from stencil
  """
  print("Sarting iteration ...")
  # Subset gap region
  jl1 = jdm - 2*dlty
  jr1 = jl1

  il1 = 0
  il2 = imid - 1
  ir1 = 2*imid - il1 - 1
  ir2 = 2*imid - il2 - 1

  AL = A2d[jl1:, il1:il2+1]
  AR = A2d[jr1:, ir2:ir1+1]

  AL_flip = np.flipud(np.fliplr(AL))
  AW = np.concatenate((AR, AL_flip), axis=0)

  AMsk = np.isnan(AW)
  #AW = mmom6.fill_land3d(AW)
  #AW = np.where(AMsk, 0.0, AW) # land points are not used 

  mm, nn = AW.shape
  mm_half = mm // 2

  # Define gap zone avoiding land:
  gap_mask = np.zeros_like(AW, dtype=bool)
  gap_mask[mm_half - nfix: mm_half + nfix, :] = True

  # do not update land or outside-gap region
  update_mask = gap_mask & (~AMsk)

  #II, JJ = np.meshgrid(np.arange(nn), np.arange(mm))

  iter = 0
  eps0 = 1.e6
  while eps0 > eps_stop and iter < Niter:
    iter += 1
    AW_old = AW.copy()  
 
    for ii in range(1, nn-1):
      for jj in range(1, mm-1):

        if not update_mask[jj,ii]:
          continue

        # For Laplace, avoid "holes" in the domain:
        # otherwise, the solution will prolagate to / from land
        w_sum = 0.
        n_sum = 0.

        if not AMsk[jj+1, ii]:
          w_sum += AW[jj+1, ii]
          n_sum += 1
        if not AMsk[jj-1, ii]:
          w_sum += AW[jj-1, ii]
          n_sum += 1
        if not AMsk[jj, ii+1]:
          w_sum += AW[jj, ii+1]
          n_sum += 1
        if not AMsk[jj, ii-1]:
          w_sum += AW[jj, ii-1]
          n_sum += 1

        if n_sum > 0:
          AW[jj,ii] = w_sum / n_sum

    # Calc. change only over the valid region
    diff = np.abs(AW - AW_old)
    eps0 = np.max(diff[update_mask])

    print(f"iter {iter:4d} | eps = {eps0:.6e}")

  # Unpack and insert back to A2d:
  #AW[AMsk] = np.nan
  AR_smth = AW[:mm_half, :]
  AL_smth = np.fliplr(np.flipud(AW[mm_half:, :]))
  A2d[jl1:, il1:il2+1] = AL_smth
  A2d[jr1:, ir2:ir1+1] = AR_smth

  return A2d
         


# Set of points along the seam line in the right half of the grid:
# Left half flips up-down and goes "on top" of the right half in polar projection
imid = idm // 2

print("Start interpolation ...")

wt_bary = None
prev_Xp = None

#dlty = 11     # N of points from the N bndry seam - defines degree of poly  = 2*(dlty - nfix)
#nfix = 6      # N of points to fix from the seam
dlty = 9     # N of points from the N bndry seam - defines degree of poly  = 2*(dlty - nfix)
nfix = 5      # N of points to fix from the seam
jj0 = jdm - 1

for il in range(imid):
  if HH[-1,il] >= 0:
    continue

  # Symmetric column index assumed
  ir = 2*imid - il - 1
  assert 0 <= ir < idm, f"Invalid symmetric index ir={ir}"

  # section idices in the Left and righ halves of the grid 
  JL = np.arange(jj0-dlty, jdm-nfix)
  IL = np.zeros_like(JL) + il

  JR = np.flipud(JL)
  IR = np.zeros_like(JR) + ir

  # Distances along the sections:
  DL = JL - jdm
  DR = jdm - JR - 1
  Xp = np.concatenate([DL, DR])

  # Interpolation target points:
  XI = np.arange(JL[-1]-jdm+1, jdm-JR[0]-1)  
  assert XI.size > 0, f"XI is empty check JL={JL} and JR={JR}"

  # Barycentr. weights compute once if Xp does not change
  if wt_bary is None:
    wt_bary = mintrp.bary_weights(Xp)
    prev_Xp = Xp.copy()
  else:
    assert np.array_equal(prev_Xp, Xp), "Xp changed between iterations"

  # Interpolate
  Yhi = interp_fld(hice,  JL, IL, JR, IR, idm, jdm, Xp, XI, wt_bary) 
  Yai = interp_fld(aice,  JL, IL, JR, IR, idm, jdm, Xp, XI, wt_bary) 
  Yhs = interp_fld(hsnow, JL, IL, JR, IR, idm, jdm, Xp, XI, wt_bary) 

  # Ensure enough values returned
  assert len(Yhi) == 2*nfix, "Not enough interpolated values (hi)"
  assert len(Yai) == 2*nfix, "Not enough interpolated values (aice)"
  assert len(Yhs) == 2*nfix, "Not enough interpolated values (hs)"

  
  hice_new[-nfix:,il] = Yhi[:nfix]
  hice_new[-nfix:,ir] = np.flipud(Yhi)[:nfix]

  aice_new[-nfix:,il] = Yai[:nfix]
  aice_new[-nfix:,ir] = np.flipud(Yai)[:nfix]

  hsnow_new[-nfix:,il] = Yhs[:nfix]
  hsnow_new[-nfix:,ir] = np.flipud(Yhs)[:nfix]

# Now solve Laplacian grad^2(A) = 0 to smooth interpolated values
print("Start Laplace smoothing")
Niter = 50
hice_new  = solve_laplace(hice_new, imid, dlty, nfix, Niter, eps_stop=0.05)
hsnow_new = solve_laplace(hsnow_new, imid, dlty, nfix, Niter, eps_stop=0.005)
aice_new  = solve_laplace(aice_new, imid, dlty, nfix, Niter, eps_stop=0.02)

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

if fsave:
  base = os.path.splitext(flice)[0]
  flrst_out = f"{base}.seam_crct.nc"

  pth_out = pthice
  dflrst_out = os.path.join(pth_out,flrst_out)
  print(f"Saving RTOFS CICE4  --> {dflrst_out}")
  ds_out.to_netcdf(dflrst_out, encoding={var: {'_FillValue': None} for var in ds_out.data_vars}, format='NETCDF3_64BIT')
  ds_out.close()

else:
  print(f"   Corrected RTOFS CICE4 is not saved ")


check_plt = False
if check_plt:
  plt.ion()

  il = 1200
  jj0 = jdm - 1
  JLL = np.arange(jj0-dlty, jdm).astype(int)
  ILL = np.zeros_like(JLL) + il

  ir = 2*imid - il - 1
  JRR = np.flipud(JLL)
  IRR = np.zeros_like(JRR) + ir

  # Distances along the sections:
  DLL = JLL - jdm
  DRR = jdm - JRR - 1

  IX = np.concatenate((DLL,DRR))

  fig1 = plt.figure(1,figsize=(9,9))
  plt.clf()

  AA = hice
  AN = hice_new.squeeze()
  Yold = np.concatenate((AA[JLL,ILL], AA[JRR,IRR]))
  Ynew = np.concatenate((AN[JLL,ILL], AN[JRR,IRR]))

  ax1 = plt.axes([0.08, 0.7, 0.8, 0.25])
  sttl = f"hice, iL={il} - iR={ir}, baryc. Lagr. interp. polynom."
  ax1.plot(IX, Yold, '.-')
  ax1.plot(IX, Ynew, '.-')
  ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
  ax1.set_xlim([IX[0]-0.5, IX[-1]+0.5])
  ax1.set_title(sttl) 

  AA = aice
  AN = aice_new.squeeze()
  Yold = np.concatenate((AA[JLL,ILL], AA[JRR,IRR]))
  Ynew = np.concatenate((AN[JLL,ILL], AN[JRR,IRR]))

  ax2 = plt.axes([0.08, 0.39, 0.8, 0.25])
  sttl = f"aice, iL={il} - iR={ir}, baryc. Lagr. interp. polynom."
  ax2.plot(IX, Yold, '.-')
  ax2.plot(IX, Ynew, '.-')
  ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
  ax2.set_xlim([IX[0]-0.5, IX[-1]+0.5])
  ax2.set_title(sttl) 

  AA = hsnow
  AN = hsnow_new.squeeze()
  Yold = np.concatenate((AA[JLL,ILL], AA[JRR,IRR]))
  Ynew = np.concatenate((AN[JLL,ILL], AN[JRR,IRR]))

  ax3 = plt.axes([0.08, 0.08, 0.8, 0.25])
  sttl = f"hsnow, iL={il} - iR={ir}, baryc. Lagr. interp. polynom."
  ax3.plot(IX, Yold, '.-')
  ax3.plot(IX, Ynew, '.-')
  ax3.grid(True, which='both', linestyle='--', linewidth=0.5)
  ax3.set_xlim([IX[0]-0.5, IX[-1]+0.5])
  ax3.set_title(sttl) 

  btx = 'correct_RTOFS_NPole_seam.py'
  bottom_text(btx, pos=[0.1, 0.04])


  #   Show section on polar projection:
  fig1 = plt.figure(1,figsize=(9,9))

  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
  rmax = 3.
  clrmp.set_bad(color=[0.1, 0.1, 0.1])
  cntr_clr = [0.9,0.,1]
   
  regn = 'north' 
  if regn == 'south':
    m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
    parallels = np.arange(-80,-10,10.)
    meridians = np.arange(-360,359.,45.)
  elif regn == 'north':
    m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
    parallels = np.arange(50, 90, 5)
    meridians = np.arange(-360, 359., 45.)

  xh, yh = m(LON, LAT) # GFS coords

  AA = hice.squeeze()     # original fields
  AA = hice_new.squeeze() # corrected fields

  plt.clf()
  ax1 = plt.axes([0.08, 0.1, 0.8, 0.8])
  #m.drawcoastlines()
  m.drawparallels(parallels,labels=[0,0,0,0],fontsize=10)
  m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10)

  img = ax1.pcolormesh(xh, yh, AA, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')

  # Plot section
  fsect = True
  if fsect:
    ax1.plot(xh[JRR,IRR], yh[JRR,IRR], '.', color=[0.5,0.5,0.5])
    ax1.plot(xh[JLL,ILL], yh[JLL,ILL], '.', color=[0.,0.9,0.8])

  sttl = f"ithkn, iL={il} - iR={ir}"
  ax1.set_title(sttl)

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


  bottom_text(btx, pos=[0.2, 0.01])



