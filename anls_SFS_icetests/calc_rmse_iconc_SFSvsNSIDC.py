"""
  RMSE of ice conc. btw SFS experiments and NSIDC NRT fields

  NSIDC fields from 
  https://noaadata.apps.nsidc.org/NOAA/G02202_V6/north/daily/2025/

  example:
  plot control run and "global" expt with inserted snow, show Arctic ocean (north)
  show persistence scores for both control and expt 21 initial states:
  run calc_rmse_iconc_NSIDC.py --regn north --enmb 1 21 --prst 1 21
  
  Show stat for twn expts (ai+hs+hsU) with old and new ("global") versions of the
  insertion codes for S. Ocean and control + persistence
  run calc_rmse_iconc_NSIDC.py --regn south --enmb 1 24 11 --prst 1 24    

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib  
import xarray
import datetime
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
init_hr = 0
expt_name = 'sfs_C192mx025_cice_test'

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere: north or south", type=str, required=True)
parser.add_argument("--init", help=f"Start of the f/cast YYYYMMDD, default=f{init_date}", choices=[20240701], type=int)
parser.add_argument("--dend", help="End date to plot YYYYMMDD or provide --ndays", type=int)
parser.add_argument("--ndays", help=f"Optional: N days to show from init, will override dend", type=int)
parser.add_argument("--pinit", help="Show initial state if exists: 1=yes (default), 0=no", default=1, type=int)
parser.add_argument(
    "--enmb",
    help="List of experiment numbers (e.g., 1 3 9 12)",
    type=int,
    nargs="+",             # <-- allows one or more integers and will generate a list
    required=True
)
parser.add_argument(
    "--prst", 
    help=f"List of experiments to show persitance, default=none",
    type=int,
    nargs="+"
)

args = parser.parse_args()
regn      = args.regn if args.regn else None
init_date = args.init if args.init else init_date  
end_date  = args.init if args.dend else None
ndays     = args.ndays if args.ndays else None
pinit     = args.pinit
ENMBS     = args.enmb if args.enmb else None
PRST      = args.prst if args.prst else []
plt_init = pinit == 1  # show RMSE for init state if init. state file exists and saved by CICE6

if end_date is None and ndays is None:
  raise RuntimeError("Both end_day and ndays are None, one of them has to be provided")

# Error in NSIDC ice concentration fields:
if regn == 'north':
  NSIDC_err = mtime.datenum_v2([[2024,7,12]], ref_day0=False)
else:
  NSIDC_err = np.array([])
 
# Dates:
# Assumed init hour = 0
if init_hr > 0:
  raise RuntimeError(f"Assumed init_hr 0, given init_hr={init_hr}, need to change code logic")
#dnmbS = mtime.rdate2datenum(init_date*100 + init_hr)  # init. day nmb
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


pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]

RMsk = np.where(HH>=0, 0, 1)
if regn == 'south':
  RMsk = np.where(hlat > -60., 0, RMsk)
elif regn == 'north':
  RMsk = np.where(hlat < 50, 0, RMsk)

def rmse2d(AA, AI):
  """
    RMSE of 2d fields
  """
  sqerr = (AA-AI)**2
 # Jice, Iice = np.where(~np.isnan(AA) & ~np.isnan(AI))
  Jice, Iice = np.where(((AA > 1.e-3) | (AI > 1.e-3)) & ~np.isnan(AA) & ~np.isnan(AI))

  Nice = len(Jice)
  assert(Nice>0), f"No ice grid found at {YR}/{MM:02d}/{DD:02d}"
  rmse_mo = np.sqrt(1./float(Nice)*np.sum(sqerr[Jice,Iice]))

  return rmse_mo


# Create an array of day numbers with 0hr = init cond, 12 hr - daily means
# Assumed: runs start at 0 hr, if not - may need to change the logic 
# for finding the ic fields 
#dnmbS = int(mtime.datenum_v2([YR,MMS,DDS], ref_day0=False))
#dnmbE = int(mtime.datenum_v2([YR,MME,DDE], ref_day0=False))
RECS = [dnmbS] + [x + 0.5 for x in range(dnmbS, dnmbE + 1)]
RECS = np.array(RECS)

nprst  = len(PRST)
nexpts = len(ENMBS)
nrecs  = RECS.shape[0]
RMSE = np.zeros((nrecs,nexpts))
if nprst > 0:
  RMSEp = np.zeros((nrecs, nprst))
Aprst = None

iens = -1
iprst = -1
for enmb in ENMBS:
  iens += 1
  irec = -1
  pthout_cice = pths_ufs[node_nm]["MOM6"]["pthcice"].format(expt_name=expt_name, enmb=enmb)

  for nn in range(nrecs):
    dnmb = RECS[nn]
    YR,MM,DD,hr = mtime.datevec(dnmb, round_hrs=True)[:4] 

    track_prst = (nprst > 0) and np.isin(enmb, PRST)
    if hr == 0:
      # Initial state
      nsec0 = 0 
      flcice = f"iceh_ic.{YR}-{MM:02d}-{DD:02d}-{nsec0:05d}.nc"
      dflcice = os.path.join(pthout_cice,flcice)

      if not os.path.isfile(dflcice) or not plt_init:
        print(f"Initial state file is missing or not requested plt_init, proceed without it ...")
        irec += 1
        RMSE[irec, iens] = np.nan
        continue
    else:
      print(f"Processing {YR}/{MM}/{DD} expt{enmb:02d}...")
      flcice = f"iceh.{YR}-{MM:02d}-{DD:02d}.nc"
      dflcice = os.path.join(pthout_cice,flcice)

    print(f"Processing {dflcice}")
    with xarray.open_dataset(dflcice) as dcice:
      AA = dcice['aice_d'].data.squeeze()

    if int(dnmb) in NSIDC_err:
      # Error ice conc fields
      rmse_mo = np.nan
    else:
      # Interpolated NSIDC obs fields:
      pthnsidc = os.path.join(pthdata,f"NRT_NOAA_NSIDC_seaconc/{YR}")
      fliceout = f'NSIDC_iconc_interp_mesh025_{jdm}x{idm}_{YR}{MM:02d}_{regn}.nc'
      dfliceout = os.path.join(pthnsidc,fliceout)
      print(f'Loading interpolated ice conc {dfliceout}')
      with xarray.open_dataset(dfliceout) as dsint:
        AI = dsint['ice_conc'].isel(time=DD-1).squeeze()

      AA = np.where(RMsk == 0, np.nan, AA)
      AI = np.where(RMsk == 0, np.nan, AI)
      rmse_mo = rmse2d(AA,AI)

    irec += 1
    if track_prst:
      if Aprst is  None:
        Aprst = AA.copy()
        iprst += 1
      rmse_prst = rmse2d(Aprst, AI)
      print(f"RMSE={rmse_mo:.3f}  RMSE_prst={rmse_prst:.3f}")
      RMSEp[irec,iprst] = rmse_prst
    else:
      print(f"RMSE={rmse_mo:.3f}")

    RMSE[irec,iens] = rmse_mo

  Aprst = None
  
# Line colors:
import mod_gfs_cice_anls as mgfscice
importlib.reload(mgfscice)
CLRS = mgfscice.sens_tests_colors()

print("Plotting ...")

XT = RECS - np.floor(RECS[0])
DV0 = mtime.datevec(RECS[0])
yr0, mm0, dd0 = DV0[:3]
ndays = np.floor(XT[-1] - XT[0])
start_date = datetime.datetime(yr0, mm0, dd0)  
end_date = start_date + datetime.timedelta(days=XT[-1])
day0 = start_date
xticks = []
xtick_labels = []
N = 40  # threshold

if ndays < N:
  # Daily ticks
  day0 = start_date
  while day0 <= end_date:
    xt = (day0 - start_date).days
    if 0 <= xt <= ndays:
      xticks.append(xt)
      xtick_labels.append(day0.strftime("%m/%d/%y"))
    day0 += datetime.timedelta(days=1)
else:
  # Day 1 and 15 of the month
  while day0 <= end_date:
    # 1st of month
    xt = (day0 - start_date).days
    if 0 <= xt <= XT[-1]:
      xticks.append(xt)
      xtick_labels.append(day0.strftime("%m/%d/%y"))

    # 15th of month
    mid = day0.replace(day=15)
    xt = (mid - start_date).days
    if 0 <= xt <= XT[-1]:
      xticks.append(xt)
      xtick_labels.append(mid.strftime("%m/%d/%y"))

    # next month
    if day0.month == 12:
      day0 = day0.replace(year=day0.year+1, month=1, day=1)
    else:
      day0 = day0.replace(month=day0.month+1, day=1)


yticks = np.arange(0.,1.,0.05)
sttl = f"RMSE btw iconc NSIDC and SFS GFS expts, {regn}\n"
sttl = sttl + f"{YR}/{MMS:02d}/{DDS:02d}-{YR}/{MME:02d}/{DDE:02d}"

plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.4, 0.8, 0.5])

LNS = []
for iens in range(nexpts):
  enmb = ENMBS[iens]
  rmse0 = RMSE[:,iens]
  clr0  = CLRS[iens,:]
  #line_lbl  = f"expt{enmb:02d}"
  line_lbl = mgfscice.sfs_tests_info(enmb)
  ln1, = ax1.plot(XT,rmse0, '.-', linewidth=2, color=clr0, label=line_lbl)
  LNS.append(ln1)

if nprst > 0:
  for iprst in range(nprst):
    enmb = PRST[iprst]
    rmse0 = RMSEp[:,iprst]
    iens = ENMBS.index(enmb)
    assert iens >= 0, f"Failed to find expt nmb for {enmb}"
    clr0  = CLRS[iens,:]
    line_lbl = f"persist expt{enmb:02d}"
    ln1, = ax1.plot(XT, rmse0, '--', linewidth=2, color=clr0, label=line_lbl)
    LNS.append(ln1)
 

yl1 = 0
yl2 = np.nanmax(RMSE) * 1.05 
#yl2 = 0.5

if nprst > 0:
  ylP = np.nanmax(RMSEp) * 1.05
  yl2 = np.max([yl2, ylP])
 
ax1.set_yticks(yticks)
ax1.set_ylim(yl1, yl2)
#ax1.set_xlim(XT[0], XT[-1])
ax1.set_xticks(xticks)
ax1.set_xticklabels(xtick_labels, rotation=60, ha='right')
ax1.grid('on')
#ax1.set_xlabel('Forecast days')
ax1.set_title(sttl)

ax3 = plt.axes([0.08, 0.1, 0.6, 0.2])
lgd = plt.legend(handles=LNS, loc='lower left')
ax3.axis('off')

btx = 'calc_rmse_iconc_SFSvsNSIDC.py'
bottom_text(btx, pos=[0.1,0.05])

