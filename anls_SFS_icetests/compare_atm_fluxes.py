"""
  Compare atm. fluxes on ice / snow
  all fluxes (excpet a few) are W/m2 averaged over ice area:
  Flux = sum(cat=1,...,nc)(flux_n(cat)*aicen(cat)) / sum(aicen(cat)) - averged over day (or requested T)

  flux_ai is (supposedly) = mean of sum(cat=1,...,nc)(flux_n(cat)*aicen(cat))

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
                   
# Append custom module paths
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

init_date = 20240701
init_hr = 0    # nominal hr, actual: -6 hrs for IAU, and -3 FHROT (f/cast hr rotation)
regn = 'north'

# All fluxes are positive downward
# fswint - sh/wave absorbed in ice interior fswint_ai_d, W/m2 <-- all 0's ?? not used
# fswabs - snow/ice/ocn absorbed solar flux
# flat   - latent heat flux, W/m2
# fsens - sens. heat flux, fsens_ai_d
# fsurf - net surface heat flux (ice to atm), positive downward, excludes conductive flux, weighted by area
# flwdn - down. longwave flux, flwdn_d
# flwup - upward longwave
# flwnet - net longwave (flwdn + flwup)
# fswdn  - downward solar (sh/wave) flux, fswdn_d
# fswup  - upward solar flux
# fswnet - net sh/wave (fswdn + fswup)
# fswthru - SW through ice to ocean (positive to the ocean)
# fhocn - turbulent ocean-to-ice heat f(Tfrz)
# fbot  - net bottom heat flux at bottom (fhocn - fcond) - what actually goes to basal melt/growth

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str)
parser.add_argument("--init", help=f"init date", 
                    choices=[20240701, 20250701], 
                    default=init_date, type=int)
parser.add_argument("--ihr", help=f"init hour, default={init_hr}", type=int)
parser.add_argument("--dend", help="End date to plot YYYYMMDD or provide --ndays", type=int)
parser.add_argument("--ndays", help=f"Optional: N days to show from init, will override dend", type=int)
parser.add_argument(
    "--fld",
    help="Plot field",
    choices=[
        'fswabs', 'flat', 'fsens', 'fsurf',
        'flwdn', 'flwup', 'flwnet',
        'fswdn', 'fswup', 'fswnet',
        'fswthru', 'fhocn', 'fbot'
    ],
    type=str,
    required=True
)
parser.add_argument(
    "--enmb",
    help="List of experiment numbers (e.g., 1 3 9 12)",
    type=int,
    nargs="+",             # <-- allows one or more integers and will generate a list
    required=True
)
args = parser.parse_args()

regn      = args.regn if args.regn else regn
init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
ENMBS     = args.enmb if args.enmb else None
end_date  = args.init if args.dend else None
ndays     = args.ndays if args.ndays else None
fld_name  = args.fld

# Get date:
plot_init = False  # initial conditions

if end_date is None and ndays is None:
  raise RuntimeError("Both end_day and ndays are None, one of them has to be provided")

dnmbS = int(mtime.rdate2datenum(init_date*100))  # init. day nmb
YRS, MMS, DDS = mtime.datevec(dnmbS)[:3]

if ndays is not None:
  # ndays overrides end_date
  dnmbE = dnmbS + ndays
else:
  dnmbE = mtime.rdate2datenum(end_date*100)
YRE, MME, DDE = mtime.datevec(dnmbE)[:3]

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
    
fyaml = 'paths_sfs.yaml'
with open(fyaml) as ff:
  pths_sfs = safe_load(ff)
    
pthgrid = pths_sfs[node_nm]["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")
    
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

# Analysis region:
RMsk = (HH < 0) 
if regn == 'north':
  RMsk &= (hlat > 50)
elif regn == 'south':
  RMsk &= (hlat < -50)

# Daily mean output time is assumed except for the initial field
if plot_init:
  RECS = [dnmbS] + [x + 0.5 for x in range(dnmbS, dnmbE + 1)]
else:
  RECS = [x + 0.5 for x in range(dnmbS, dnmbE + 1)]

aice_eps = 0.95 # minimum ice conc to consider
hbins = np.array([0,0.5,1,2,100])
ncats = len(hbins)-1

RECS   = np.array(RECS)
nrecs  = RECS.shape[0]
nexpts = len(ENMBS)
FMLT = np.zeros((nexpts, nrecs, ncats))

#cnvrt_ice_mean = False  # convert to ice area mean - not needed as all fluxes should be wrt ice area
avrg = "_ice"
varnm2 = None
swap_sign = False
if fld_name == 'fswabs':
  varnm = 'fswabs_d'     # snow/ice/ocn absorbed solar flux
  strs = "Sn/ice/ocn absorbed solar, W/m2"
elif fld_name == 'flat':
  varnm = 'flat_d'    # latent heat flux (cpl), positive down, cpl - cell averaged?
  strs = "Latent heat flux (cpl), W/m2"
elif fld_name == 'fsens':
  varnm = 'fsens_d'       # sens heat flux m2/m2_ice_area
  strs = f"Sens flux, W/m2 "
  cnvrt_ice_mean = True
elif fld_name == 'fsurf':
  varnm = 'fsurf_ai_d' # positive downward, excludes conductive flux, weighted by area
  strs = f"Net surf heat flux, W/m2 "
  avrg = '_cell'
elif fld_name == 'flwdn':
  varnm = 'flwdn_d'     # down longwave, positive downward
  strs = "Down Longwave, W/m2"
elif fld_name == 'flwup':
  varnm = 'flwup_d'     # upward longwave (cpl), positive downward
  strs = "Upward Longwave, W/m2"
elif fld_name == 'flwnet':
  varnm = 'flwup_d'
  varnm2 = 'flwdn_d'
  strs = 'Net longwave, W/m2' 
elif fld_name == 'fswdn':
  varnm = 'fswdn_d'      # downward solar (sh/wave) flux, positive downward
  strs = "Down sh/wave, W/m2"
elif fld_name == 'fswup':
  varnm = 'fswup_d'      # Upward solar (sh/wave) flux, positive upward
  swap_sign = True       # positive downward, i.e. reflected sh/wave < 0
  strs = "Upward sh/wave, W/m2"
elif fld_name == 'fswnet':
  varnm = 'fswup_d'      # Upward solar (sh/wave) flux, positive up
  varnm2 = 'fswdn_d'     # Downward solar (sh/wave) flux, positive down
  swap_sign = True
  strs = "Net sh/wave, W/m2"
elif fld_name == 'fswthru':
  varnm = 'fswthru_d'  # Sh/wave through ice to ocean (cpl)
  strs = "Sh/wave through ice to ocean, W/m2"
elif fld_name == 'fhocn':
  varnm = 'fhocn_d'   # heat flux ice to ocn (cpl)
  strs = "Heat flux ice to ocn, W/m2"
elif fld_name == 'fbot':
  varnm = 'fbot_d'  # ice net heat flux = ocean_to_ice - Fcond --> basal ice melt/growth
  strs = "Basal net heat flux, W/m2"

strs = strs + avrg  # m2 of ice area or grid cell area

iens = -1
for enmb in ENMBS:
  iens += 1
  irec = -1
  pthout_cice = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/sfs_C192mx025_cice_test/expt{enmb:02d}/cice6"

  for nn in range(nrecs):
    irec += 1
    dnmb = RECS[nn]
    YR, MM, DD, hr = mtime.datevec(dnmb, round_hrs=True)[:4]

    if hr == 0:
      # Initial state
      nsec0 = 0
      flcice = f"iceh_ic.{YR}-{MM:02d}-{DD:02d}-{nsec0:05d}.nc"
      dflcice = os.path.join(pthout_cice,flcice)

      if not os.path.isfile(dflcice) or not plot_init:
        print(f"Initial state file is missing or not requested plt_init, proceed without it ...")
        RMSE[irec, iens] = np.nan
        continue
    else:
      print(f"Processing {YR}/{MM}/{DD} expt{enmb:02d}...")
      flcice = f"iceh.{YR}-{MM:02d}-{DD:02d}.nc"
      dflcice = os.path.join(pthout_cice,flcice)

    print(f"Processing {dflcice}")
    with xarray.open_dataset(dflcice) as dcice:
      AI = dcice['aice_d'].values.squeeze()
      HI = dcice['hi_d'].values.squeeze()
      AF = dcice[varnm].values.squeeze() 
      #units = dcice[varnm].attrs["units"]
      if swap_sign: 
        AF = - AF
      if varnm2 is not None:
        # net fluxes
        AF2 = dcice[varnm2].values.squeeze()
        AF += AF2

    # averaged over ice area:
    #if cnvrt_ice_mean:
    #  AF = np.divide(AF, AI, out=np.zeros_like(AF), where=AI > 0)

    AF[AI < aice_eps] = np.nan
    AF[~RMsk] = np.nan
    HI[AI < aice_eps] = np.nan
    HI[~RMsk] = np.nan

    melt_cat = np.full(ncats, np.nan)

    for ic in range(ncats):
      hmin = hbins[ic]
      hmax = hbins[ic+1]

      mask = (HI >= hmin) & (HI < hmax)

      # combine with existing valid-data mask
      mask = mask & ~np.isnan(HI) & ~np.isnan(AF)

      if np.any(mask):
        melt_cat[ic] = np.nanmean(AF[mask])

    FMLT[iens,irec,:] = melt_cat


# Line colors:
import mod_gfs_cice_anls as mgfscice
importlib.reload(mgfscice)
CLRS = mgfscice.sens_tests_colors()
XT = RECS - np.floor(RECS[0])

nbins = len(hbins) - 1
ncols = 2
nrows = nbins // ncols + nbins % ncols

if nrows <= 2:
    bottom = 0.3
elif nrows == 3:
    bottom = 0.12
else:
    bottom = 0.1

sinfo = f'SFS init {init_date}, {regn}, ' + strs + f", {varnm}," + f" aice>{aice_eps:.2f}"
#sinfo = sinfo + f'{pthout_cice}'

plt.ion()

fig1 = plt.figure(1,figsize=(9,9))
fig1.clf()  
axes = fig1.subplots(nrows=nrows, ncols=ncols)

fig1.subplots_adjust(
    left=0.05,
    right=0.99,  
    top=0.95,
    bottom=bottom,
    wspace=0.12,
    hspace=0.12
)

for ax in axes.flat:
  ax.set_box_aspect(0.6)

for ibin in range(nbins):
  hmin = hbins[ibin]
  hmax = hbins[ibin+1]

  irow = ibin // ncols
  icol = ibin % ncols
  ax1 = axes[irow, icol]

  fmlt_max = 0. 
  for iens in range(nexpts):
    enmb = ENMBS[iens]
    fmlt = FMLT[iens,:,ibin]
    clr0  = CLRS[iens,:]
    fmlt_max = np.max([np.max(fmlt), fmlt_max])
    ax1.plot(XT, fmlt, '.-', linewidth=2, color=clr0)

  ax1.grid(True, alpha=0.3)
  #ax1.set_xlabel('Forecast days')
  if hmax < 10:
    ax1.set_title(f'hice=[{hmin:.1f}, {hmax:.1f}]')
  else:
    ax1.set_title(f'hice>{hmin:.1f}')

  if 0 < fmlt_max < 0.01:
    ax1.set_ylim([0, 0.1])
   
line_lbl = mgfscice.sfs_tests_info(enmb)

ax3 = plt.axes([0.55, 0.04, 0.43, 0.05])
handles = [
    plt.Line2D([0], [0], color=CLRS[i,:], lw=2)
    for i in range(nexpts)
]

leg_labels = []
for enmb in ENMBS:
  lbl = mgfscice.sfs_tests_info(enmb)
  leg_labels.append(f"{lbl}")

ax3.legend(handles, leg_labels, loc='lower left')
ax3.axis('off')

# Overasll header:
ax4 = plt.axes([0.01, 0.95, 0.98, 0.03])
ax4.text(
    0.5, 0.5, sinfo,
    fontsize=12,
    ha='center',   # horizontal align
    va='bottom',   # vertical align
    transform=ax4.transAxes
)
ax4.axis('off')

btx = 'compare_atm_fluxes.py'
bottom_text(btx, pos=[0.01,0.015])

