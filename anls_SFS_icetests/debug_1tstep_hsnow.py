"""
  Check snow depth at a grid point

  From a debug run with saved variables in ice thickness categories
  Output - everty time step

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
import mod_misc1 as mmisc

init_date = 20250701
init_hr = 0

# 62 - debug run with every step instant. output !!! 
enmb = 62
dt = 600         # CICE6 time step, sec - check with ice_in

parser = argparse.ArgumentParser()
parser.add_argument("--init", help=f"init date", choices=[20240701, 20250701], default=init_date, type=int)
parser.add_argument("--ihr", help=f"init hour, default {init_hr}", type=int)
parser.add_argument("--fday1", help=f"forecast day to plot: 1, 2, ... ", type=int, required=True)
parser.add_argument("--nsec1", help="seconds in day1, dt=600 sec: 0, 600, 1200, ..., 85800, fday=1 and nsec=0 - init cond", 
                    type=int, required=True)
parser.add_argument("--fday2", help=f"forecast day to plot: 0, 1, 2, ..., default=fday1 ", type=int)
parser.add_argument("--nsec2", help=f"seconds in day2, 0, 600, 1200, ..., 85800, default=nsec1+{dt}", 
                    type=int)
parser.add_argument("--enmb", help="experiment number: 0, 1, 2, ...", default=enmb, type=int)
args = parser.parse_args()

init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
fday1     = args.fday1
nsec1     = args.nsec1
fday2     = args.fday2 if args.fday2 is not None else args.fday1
nsec2     = args.nsec2 if args.nsec2 is not None else args.nsec1 + dt
enmb      = args.enmb
rho_snow = 300.
tav = '1'

# Select grid point for analysis:
# Grid point analysis
# Rapid snow change:
ii0 = 303
jj0 = 1062

# Low snow change:
#ii0 = 311
#jj0 = 1070


if fday1 == 0:
  fday1 = 1

if fday2 == 0:
  fday2 = 1

assert nsec1 >= 0 and nsec1 <= 86400-dt, f"nsec1 is out of bound for seconds {nsec1}"
assert nsec2 >= 0 and nsec1 <= 86400-dt, f"nsec2 is out of bound for seconds {nsec2}"

def path_names(fday, nsec, init_date, init_hr, enmb):
  plot_init = fday == 1 and nsec == 0  # initial conditions
  
  dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
  YRI,MMI,DDI,hrI = mtime.datevec(dnmbI, round_hrs=True)[:4]

  if plot_init:
    dnmb0 = dnmbI
  else:
    dnmb0 = dnmbI + fday-1                           # day to plot

  yr0,mm0,dd0,hr0 = mtime.datevec(dnmb0, round_hrs=True)[:4]
  nsec0 = int(hr0 * 3600)
  YR,MM,DD = mtime.datevec(dnmb0)[:3]

  pthoutp = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/sfs_C192mx025_cice_test/expt{enmb:02d}/cice6"
  if plot_init:
    flinp = f"iceh_ic.{yr0}-{mm0:02d}-{dd0:02d}-{nsec0:05d}.nc"
    tav = 'd'
  else:
    flinp = f"iceh_inst.{yr0}-{mm0:02d}-{dd0:02d}-{nsec:05d}.nc"
    tav = '1'

  dflice = os.path.join(pthoutp,flinp)
  #print(f"Reading {YR}/{MM}/{DD}, SFS init {YRI}/{MMI:02d}/{DDI:02d}\n  {dflice}")

  return dflice, tav

def read_cice(fday, nsec, HH, init_date, init_hr, enmb):
  dflice, tav = path_names(fday, nsec, init_date, init_hr, enmb)
  print(f"Reading {dflice}")

  with xarray.open_dataset(dflice) as dcice:
    vsno_m2c = dcice[f'hs_{tav}'].data.squeeze()   # grid cell mean snow vol, m3/m2_cell
    aicen = dcice[f'aicen_{tav}'].values.squeeze()  # ice area in cats
    vicen = dcice[f'vicen_{tav}'].values.squeeze()  # ice vol m2_ice in cats
    vsnon = dcice[f'vsnon_{tav}'].values.squeeze()  # snow vol m2_snow (!) in cats

  vsno_m2c[HH >= 0] = np.nan

  return vsno_m2c, aicen, vicen, vsnon

def read_cicefld_gridpnt(fday, nsec, init_date, init_hr, varnm, jj0, ii0, enmb):
  dflice, tav = path_names(fday, nsec, init_date, init_hr, enmb)
  #print(f"Reading {varnm} from {dflice}")

  with xarray.open_dataset(dflice) as dcice:
    FLD = dcice[f"{varnm}_{tav}"].values.squeeze()

  if FLD.ndim == 2:
    apnt = FLD[jj0,ii0]
  elif FLD.ndim == 3:
    apnt = FLD[:,jj0,ii0]
  elif FLD.ndim == 4:
    apnt = FLD[:,:,jj0,ii0]
  else:
    raise ValueError(
        f"Unsupported dimensions for {varnm}: {FLD.shape}"
    )

  return apnt

def read_cicefld_attr(fday, nsec, init_date, init_hr, var, jj0, ii0, enmb):
  # This is for scalar output only, i.e. not for (ncat,:,:) fields
  dflice, tav = path_names(fday, nsec, init_date, init_hr, enmb)
  varnm = f"{var}_{tav}"
  #print(f"Reading {varnm} from {dflice}")

  with xarray.open_dataset(dflice) as dcice:
    if varnm not in dcice:
      #raise KeyError(f"{varnm} not found in {dflice}")
      apnt = 99999
      units = 'XXXX'
      long_name = f"{varnm} missing in output"
   
      return apnt, units, long_name
  
    FLD = dcice[varnm].values.squeeze()
    units = dcice[varnm].attrs.get('units', '')
    long_name = dcice[varnm].attrs.get('long_name', '')

  if FLD.ndim == 2:
    apnt = FLD[jj0,ii0]
  elif FLD.ndim == 3:
    apnt = FLD[:,jj0,ii0]
  else:
    raise ValueError(
        f"Unsupported dimensions for {varnm}: {FLD.shape}"
    )

  return apnt, units, long_name

def print_diagn_pnt(fday1, nsec1, fday2, nsec2, init_date, init_hr, jj0, ii0, enmb, VARS):
  lname_w = 40
  units_w = 10
  varnm_w = 15

  def fmt_val(x):
    if 0 < abs(x) < 1.e-6:
      return f"{x:12.5e}"
    else:
      return f"{x:12.7f}"

  # Print 1 pnt diagn for 2D flds (not cats)
  #sfx = '_1' # every time step
  for var in VARS:
    #varnm = f"{var}{sfx}"

    val1, _, _ = read_cicefld_attr(fday1, nsec1, init_date, init_hr, var, jj0, ii0, enmb)
    val2, units, lname = read_cicefld_attr(fday2, nsec2, init_date, init_hr, var, jj0, ii0, enmb)
    
    #print(f"Reading {varnm}: {units}  {lname}")
    # Truncate if long name is too long: 
    print(
        f"{lname[:lname_w]:<{lname_w}} "
        f"{units:<{units_w}} "
        f"{var:<{varnm_w}} :  "
        f"{fmt_val(val1)}   {fmt_val(val2)}"
    )

  return
  
syst_info = os.uname() 
machine = syst_info.nodename
  
if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea compute node:", machine)
  node_nm = "gaea"
else:
  print("Unknown machine:", machine)

fyaml = 'paths_ufs.yaml'
with open(fyaml) as ff:
  pths_ufs = safe_load(ff)

# Get MOM6 grid
pthgrid = pths_ufs[node_nm]["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')


with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

Vsno1_m2c, Aicen1, Vicen1, Vsnon1 = read_cice(fday1, nsec1, HH, init_date, init_hr, enmb)
Vsno2_m2c, Aicen2, Vicen2, Vsnon2 = read_cice(fday2, nsec2, HH, init_date, init_hr, enmb)

clrmp = mclrmps.colormap_temp()
rmin = 0.
rmax = 0.1

clrmp.set_bad(color=[0.1, 0.1, 0.1])
cntr_clr = [0.9,0.,1]

plt.ion()

print("Plotting ...")


A2d = Vsno2_m2c  # grid cell mean snow depth (vol snow / m2_grid)
#A2d = aggr_vsn2


nstep = nsec2 // dt
sttl = f"SFS expt{enmb} init {init_date} FDAY={fday2} nsec={nsec2}, tstep={nstep}\ngrid cell mean hsnow"

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.08, 0.1, 0.8, 0.8])

img = ax1.pcolormesh(A2d, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')

ax1.set_aspect('equal', adjustable='box')
ax1.set_ylim(800, 1080)
ax1.set_xlim(100, 670)
ax1.set_title(sttl)

ax1.plot(ii0, jj0, marker='o', color=[1, 0.0, 0.8])

ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax3.set_xticklabels(ax3.get_xticks())
ticklabs = clb.ax.get_xticklabels()
clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=12)
clb.ax.tick_params(direction='in', length=12)


btx = 'debug_1step_hsnow.py'
bottom_text(btx)


# Get values in the grid
hs1 = Vsno1_m2c[jj0,ii0]  # Grid cell mean snow thickness = m3(snow) / m2_cell
hs2 = Vsno2_m2c[jj0,ii0]
aicen1 = Aicen1[:,jj0,ii0]
aicen2 = Aicen2[:,jj0,ii0]
vicen1 = Vicen1[:,jj0,ii0]
vicen2 = Vicen2[:,jj0,ii0]
vsnon1 = Vsnon1[:,jj0,ii0]
vsnon2 = Vsnon2[:,jj0,ii0]


aice1 = np.sum(aicen1)
aice2 = np.sum(aicen2)
hicen1 = vicen1 / aicen1
hicen2 = vicen2 / aicen2
vsno1 = np.sum(vsnon1 * aicen1)  # m3/m2_ice
vsno2 = np.sum(vsnon2 * aicen2)

VARS = ['aice', 'hi', 'hs', 'snowfrac', 'Tair', 'Tsfc', 'sice', 
        'fswdn', 'fswup','flwdn', 'flwup', 'snow', 'rain', 'frzmlt', 'albsni',
        'albsno', 'albice', 'albpnd', 'flat', 'fsens', 'fsurf_ai', 'congel',
        'frazil', 'snoice', 'dsnow', 'melts', 'meltt',
        'meltb', 'meltl', 'fbot', 'fhocn', 'fswthru', 'apond',
        'hpond', 'ipond',  'fsloss','rsnw', 'smassliq', 'meltsliq']

print(f" ------------------------")
print(f"Grid pnt: j={jj0} i={ii0}\n HH={HH[jj0,ii0]:.2f}, lon={hlon[jj0,ii0]:.2f} lat={hlat[jj0,ii0]:.2f}")
print_diagn_pnt(fday1, nsec1, fday2, nsec2, init_date, init_hr, jj0, ii0, enmb, VARS)

# 3D fields with ice cats.
str = f"vsnon j={jj0} i={ii0}:"
mmisc.print_cols(vsnon1, vsnon2, str=str)

str = f"aicen j={jj0} i={ii0}:"
mmisc.print_cols(aicen1, aicen2, str=str)

str = f"vicen j={jj0} i={ii0}:"
mmisc.print_cols(vicen1, vicen2, str=str)

str = f"hice j={jj0} i={ii0}:"
mmisc.print_cols(hicen1, hicen2, str=str)



