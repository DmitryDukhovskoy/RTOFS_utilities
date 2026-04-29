"""
  Plot ice bottom melt 
  from datm experiment

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
import mod_gfs_cice_anls as mgfscice
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

YR = 2025
MMS = 1 
DDS = 3
MME = 1
DDE = 16
regn = 'south'

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str)
parser.add_argument("--yr", help=f"year of model run, default={YR}", type=int)
parser.add_argument("--ms", help=f"start month of data to plot, default={MMS}", type=int)
parser.add_argument("--me", help=f"end month of NSIDCS data to plot, default={MME}", type=int)
parser.add_argument("--ds", help=f"Start: day in the start month, default={DDS}", type=int)
parser.add_argument("--de", help=f"End: day in the end month, default={DDE}", type=int)
parser.add_argument(
    "--enmb",
    help="List of experiment numbers (e.g., 1 3 9 12)",
    type=int,
    nargs="+",             # <-- allows one or more integers and will generate a list
    required=True
)
args = parser.parse_args()
  
regn  = args.regn if args.regn else regn
YR    = args.yr if args.yr else YR
MMS   = args.ms if args.ms else MMS
MME   = args.me if args.me else MMS
DDS   = args.ds if args.ds else DDS
DDE   = args.de if args.de else DDE
ENMBS = args.enmb if args.enmb else None
  
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

pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthnsidc = os.path.join(pthdata,f"NRT_NOAA_NSIDC_seaconc/{YR}")
if regn == 'south':
  RMsk = np.where(HH>=0, 0, 1)
  RMsk = np.where(hlat > -60., 0, RMsk)

# Create an array of day numbers with 0hr = init cond, 12 hr - daily means
dnmbS = int(mtime.datenum([YR,MMS,DDS]))
dnmbE = int(mtime.datenum([YR,MME,DDE]))
RECS = [x + 0.5 for x in range(dnmbS, dnmbE + 1)]
RECS = np.array(RECS)

nexpts = len(ENMBS)
nrecs  = RECS.shape[0]
MELTB = np.zeros((nrecs,nexpts))
iens = -1
for enmb in ENMBS:
  iens += 1
  irec = -1
  pthout_cice = pths_ufs[node_nm]["MOM6"]["pthcice"].format(enmb=enmb)

  for nn in range(nrecs):
    dnmb = RECS[nn]
    YR,MM,DD,hr = mtime.datevec(dnmb, round_hrs=True)[:4] 

    print(f"Processing {YR}/{MM}/{DD} expt{enmb:02d}...")
    flcice = f"iceh.{YR}-{MM:02d}-{DD:02d}.nc"
    dflcice = os.path.join(pthout_cice,flcice)

    print(f"Processing {dflcice}")
    with xarray.open_dataset(dflcice) as dcice:
      AA = dcice['meltb_d'].data.squeeze()      # basal ice melt cm/day
      aice = dcice['aice_d'].data.squeeze()

    AA = np.where(RMsk == 0, np.nan, AA)
    AA = np.where(aice < 0.15, np.nan, AA)   # ignore near ice edge regions

    melt_mn = np.nanmean(AA)
    
    print(f"melt={melt_mn:.4f} cm/day")
    irec += 1
    MELTB[irec,iens] = melt_mn


# Line colors:
CLRS = np.array([[0., 0.2, 0.9],
                 [0.7, 0., 1],
                 [0., 0.8, 0.3],
                 [0., 0.8, 1],
                 [0.9, 0.4, 0],
                 [1., 0., 0],
                 [0.5, 0.3, 0],
                 [0.5, 0.5, 0.5],
                 [0.7, 0.45, 0.9]])


print("Plotting ...")

XT = RECS - np.floor(RECS[0])
xticks = np.arange(np.floor(XT[0]),np.ceil(XT[-1]))
sttl = f"Ice basal melt (cm/day) in datmUFS expts, {YR}/{MMS:02d}/{DDS:02d}-{YR}/{MME:02d}/{DDE:02d}"

plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.4, 0.8, 0.5])

LNS = []
for iens in range(nexpts):
  enmb = ENMBS[iens]
  mltb0 = MELTB[:,iens]
  clr0  = CLRS[iens,:]
  #line_lbl  = f"expt{enmb:02d}"
  line_lbl = mgfscice.sens_tests_info(enmb)
  ln1, = ax1.plot(XT,mltb0, 'o-', linewidth=2, color=clr0, label=line_lbl)
  LNS.append(ln1)

ax1.set_xticks(xticks)
ax1.grid('on')
ax1.set_xlabel('Forecast days')
ax1.set_title(sttl)

ax3 = plt.axes([0.1, 0.15, 0.6, 0.2])
lgd = plt.legend(handles=LNS, loc='upper left')
ax3.axis('off')

btx = 'icebottom_melt_datmUFS.py'
bottom_text(btx, pos=[0.1,0.1])

