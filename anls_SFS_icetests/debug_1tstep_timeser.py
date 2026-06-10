"""
  Check snow depth at a grid point

  From a debug run with saved variables in ice thickness categories
  Output - everty time step

FIelds available for diagnostics:
rate of snow loss to leads (liquid)      kg/m^2/s   fsloss          :     0.0000000    2.16454e-21
ice area  (aggregate)                    1          aice            :  
grid cell mean ice thickness             m          hi              :     1.3699065      1.3729615
grid cell mean snow thickness            m          hs              :     0.0360882      0.0331190
grid cell mean snow fraction             1          snowfrac        :     1.0000000      0.9999992
air temperature                          C          Tair            :   253.0000000      1.4122521
snow/ice surface temperature             C          Tsfc            :    -1.0081611    -11.4902821
bulk ice salinity                        ppt        sice            :     4.2866573      4.3212566
down solar flux                          W/m^2      fswdn           :     0.0000000      0.0000000
upward solar flux                        W/m^2      fswup           :     0.0000000      0.0000000
down longwave flux                       W/m^2      flwdn           :   180.0000000      0.0000000
upward longwave flux (cpl)               W/m^2      flwup           :     0.0000000   -252.4975128
snowfall rate (cpl)                      cm/day     snow            :     0.0000000      0.0000000
rainfall rate (cpl)                      cm/day     rain            :     0.0000000      0.0000000
freeze/melt potential                    W/m^2      frzmlt          :     0.0000000      0.0000000
snow/ice broad band albedo               %          albsni          :     0.8585612      0.0000000
snow albedo                              %          albsno          :     0.8585612      0.0000000
bare ice albedo                          %          albice          :     0.0000000      0.0000000
melt pond albedo                         %          albpnd          :     0.0000000      0.0000000
latent heat flux (cpl)                   W/m^2      flat            :     0.0000000     35.4926414
sensible heat flux (cpl)                 W/m^2      fsens           :     0.0000000     84.5375366
net surface heat flux                    W/m^2      fsurf_ai        :     0.0000000   -132.4673309
congelation ice growth                   cm/day     congel          :     0.0000000      0.0012629
frazil ice growth                        cm/day     frazil          :     0.0000000      0.0000000
snow-ice formation                       cm/day     snoice          :     0.0000000     43.1153679
snow formation                           cm/day     dsnow           :     0.0000000      0.3269257
top snow melt                            cm/day     melts           :     0.0000000    2.36939e-16
top ice melt                             cm/day     meltt           :     0.0000000    6.51258e-33
basal ice melt                           cm/day     meltb           :     0.0000000      0.5831583
lateral ice melt                         cm/day     meltl           :     0.0000000      0.0000000
heat flux ice to ocean (fbot)            W/m^2      fbot            :     0.0000000      0.0000000
heat flux ice to ocn (cpl)               W/m^2      fhocn           :     0.0000000     24.1388054
SW thru ice to ocean (cpl)               W/m^2      fswthru         :     0.0000000      0.0000000
melt pond fraction of sea ice            1          apond           :     0.0000000      0.0005333
mean melt pond depth over sea ice        m          hpond           :     0.0000000    2.10736e-22
mean pond ice thickness over sea ice     m          ipond           :     0.0000000    8.62186e-17
rate of snow loss to leads (liquid)      kg/m^2/s   fsloss          :     0.0000000    2.16454e-21
average snow grain radius                10^-6 m    rsnw            :   100.0000000    100.8501358
liquid mass per unit area in snow        kg/m^2     smassliq        :  99999.0000000    5.01951e-19
snow liquid contribution to meltponds    kg/m^2/s   meltsliq        :  99999.0000000      0.0000000

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
# 
# Variables in CICE6:
#
parser = argparse.ArgumentParser()
parser.add_argument("--init", help=f"init date", choices=[20240701, 20250701], default=init_date, type=int)
parser.add_argument("--ihr", help=f"init hour, default {init_hr}", type=int)
parser.add_argument("--fday1", help=f"Start forecast day to plot: 1, 2, ... ", type=int, required=True)
parser.add_argument("--nsec1", help="seconds in day1, dt=600 sec: 0, 600, 1200, ..., 85800, fday=1 and nsec=0 - init cond", 
                    type=int, required=True)
parser.add_argument("--fday2", help=f"End forecast day to plot: 0, 1, 2, ..., default=fday1 ", type=int)
parser.add_argument("--nsec2", help="seconds in day2, 0, 600, 1200, ..., 85800", type=int, required=True)
parser.add_argument("--enmb", help=f"experiment number: 0, 1, 2, ..., default={enmb}", 
                    default=enmb, type=int)
parser.add_argument("--varnm", help="variable name to plot", required=True, type=str)
args = parser.parse_args()

init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
fday1     = args.fday1
nsec1     = args.nsec1
fday2     = args.fday2 if args.fday2 is not None else args.fday1
nsec2     = args.nsec2
enmb      = args.enmb
varnm0    = args.varnm
rho_snow = 300.
tav = '1'

if fday1 == 0:
  fday1 = 1

if fday2 == 0:
  fday2 = 1

assert nsec1 >= 0 and nsec1 <= 86400-dt, f"nsec1 is out of bound for seconds {nsec1}"
assert nsec2 >= 0 and nsec1 <= 86400-dt, f"nsec2 is out of bound for seconds {nsec2}"

# Specify Grid point for analysis
# Rapid snow change:
IJp = [
       (303, 1062),     # rapid snow change point
       (311, 1070)      # slow snow change
      ]

IJp = np.array(IJp)

# Output time stamps:
RECT = []
nsec_end_day = int(86400-dt)
for iday in range(fday1, fday2+1):
  for nsec in range(0, nsec_end_day + 1, dt):
    if iday == fday1 and nsec < nsec1:
      continue
    elif iday == fday2 and nsec > nsec2:
      break

    RECT.append((iday,nsec))

RECT = np.array(RECT)


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

  if plot_init:
    flinp = f"iceh_ic.{yr0}-{mm0:02d}-{dd0:02d}-{nsec0:05d}.nc"
  #  tav = 'd'
  else:
    flinp = f"iceh_inst.{yr0}-{mm0:02d}-{dd0:02d}-{nsec:05d}.nc"
  #  tav = '1'

  pthoutp = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/sfs_C192mx025_cice_test/expt{enmb:02d}/cice6"
  dflice = os.path.join(pthoutp,flinp)
  #print(f"Reading {YR}/{MM}/{DD}, SFS init {YRI}/{MMI:02d}/{DDI:02d}\n  {dflice}")

  # Find extensions, assuming var name = varnm0_? (d, h, 1, ...)
  with xarray.open_dataset(dflice) as ds:
    varname = next(v for v in ds.variables
                   if v.startswith(f"{varnm0}_"))
    tav = varname.split("_", 1)[1]

  return dflice, tav

def read_2Dcice(fday, nsec, HH, init_date, init_hr, enmb, var):
  dflice, tav = path_names(fday, nsec, init_date, init_hr, enmb)
  print(f"Reading {dflice}")

  with xarray.open_dataset(dflice) as dcice:
    A2d = dcice[f'{var}_{tav}'].values.squeeze()  

  A2d[HH >= 0] = np.nan

  return A2d

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



A2d = read_2Dcice(fday2, nsec2, HH, init_date, init_hr, enmb, varnm0)

# Derive time series:
npp = IJp.shape[0]
nrec = len(RECT)
VALS = np.zeros((nrec,npp)) * np.nan
for ipp in range(npp):
  ii0, jj0 = IJp[ipp, :]

  print(f"Processing ii0={ii0}, jj0={jj0}")
  for irec in range(nrec):
    fday, nsec = RECT[irec,:]

    val, units, lname = read_cicefld_attr(fday, nsec, init_date, init_hr, varnm0, jj0, ii0, enmb)
    VALS[irec, ipp] = val
  
    
# Show points:
clrmp = mclrmps.colormap_temp()
rmin = np.nanmin(A2d)
rmax = np.nanmax(A2d)
if varnm0 == 'hs':
  rmax = 0.1
  rmin = 0. 
elif varnm0 == 'snowfrac':
  rmax = 1
  rmin = 0.


clrmp.set_bad(color=[0.1, 0.1, 0.1])
cntr_clr = [0.9,0.,1]

plt.ion()

print("Plotting ...")


nstep = nsec2 // dt
sttl = f"SFS expt{enmb} init {init_date} FDAY={fday2} nsec={nsec2}, tstep={nstep}\n{varnm0}"

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.08, 0.55, 0.8, 0.4])

img = ax1.pcolormesh(A2d, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')

ax1.set_aspect('equal', adjustable='box')
ax1.set_ylim(850, 1080)
ax1.set_xlim(100, 670)
ax1.set_title(sttl)

ax1.scatter(
    IJp[:,0],
    IJp[:,1],
    s=60,
    facecolor=[1, 0.0, 0.8],
    edgecolor='k',
    linewidth=0.5
)

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
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=12)
clb.ax.tick_params(direction='in', length=12)

# Plot time series:
# Line colors:
#import mod_gfs_cice_anls as mgfscice
#importlib.reload(mgfscice)
#CLRS = mgfscice.sens_tests_colors()
# No more than 10 lines:
BASE_COLORS = [
    'tab:blue', 'tab:orange', 'tab:green', 'tab:red',
    'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray',
    'tab:olive', 'tab:cyan'
]

CLRS = BASE_COLORS[:npp]

ax3 = fig1.add_axes([0.08, 0.1, 0.9, 0.4])
LBLS = []
for ipp in range(npp):
  clr0 = CLRS[ipp]
  ii0, jj0 = IJp[ipp,:]
  ln1, = ax3.plot(VALS[:,ipp], '-', linewidth=2, color=clr0, label=f"{ipp+1}: i={ii0} j={jj0}")
  LBLS.append(ln1)

#ax3.set_xlim([0.15, 1.05])
ax3.grid('on')

sttl3 = f"{varnm0} ({units}) {lname}" 
ax3.set_title(sttl3)
ax3.set_xlabel('CICE6 Time steps')

ax3.legend(
    handles=LBLS,
    loc='lower right',
    bbox_to_anchor=(1.01, -0.25),
    fontsize=9
)

btx = 'debug_1step_timeser.py'
bottom_text(btx, pos=[0.02, 0.02])




