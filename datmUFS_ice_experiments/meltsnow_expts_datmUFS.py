"""
  Timeser of surface snow melt rates
  averaged over ice 

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib  
import xarray
from copy import copy
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
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

YR = 2025
MM = 1 
DD = 3
DDE = 16

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere: north or south", type=str, required=True)
parser.add_argument("--yr", help=f"year of model run, default={YR}", type=int)
parser.add_argument("--ms", help=f"start month of data to plot, default={MM}", type=int)
parser.add_argument("--me", help=f"end month of NSIDCS data to plot, default={MM}", type=int)
parser.add_argument("--ds", help=f"Start: day in the start month to plot, default={DD}", type=int)
parser.add_argument("--de", help=f"End: day in the end month to plot, default={DDE}", type=int)
parser.add_argument(
    "--enmb",
    help="List of experiment numbers (e.g., 1 3 9 12)",
    type=int,
    nargs="+",             # <-- allows one or more integers and will generate a list
    required=True
)
args = parser.parse_args()
  
regn  = args.regn if args.regn else None
YR    = args.yr if args.yr else YR
MMS   = args.ms if args.ms else MM
MME   = args.me if args.me else MMS
DDS   = args.ds if args.ds else DD
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

DX, DY = mmom6.dx_dy(hlon, hlat)
Acell = DX*DY


with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

jdm, idm = HH.shape

pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthnsidc = os.path.join(pthdata,f"NRT_NOAA_NSIDC_seaconc/{YR}")

RMsk = np.where(HH>=0, 0, 1)
if regn == 'south':
  RMsk = np.where(hlat > -60., 0, RMsk)
else:
  RMsk = np.where(hlat < 50, 0, RMsk)

def spavrg_field(AA, HS, Acell, hsn_min=0.05, alfa=25):
  """
    Spatial average of 2d fields
    weighted by cell area
 
    hsn_min - min snow depth where melt is checked
              if hsn_min is too low, this artificially 
              reduces the average snow melt (dh/dt < hsnow, cannot melt more)
  """
  Jice, Iice = np.where( (HS > 0.05) & ~np.isnan(AA) )
  Nice = len(Jice)
  assert(Nice>0), f"No ice grid found at {YR}/{MM:02d}/{DD:02d}"
  area_tot = np.sum(Acell[Jice,Iice])

  AVRG = np.sum(AA[Jice,Iice]*Acell[Jice,Iice]) / area_tot
  Uprc = np.percentile(AA[Jice,Iice], (100-alfa))
  Lprc = np.percentile(AA[Jice,Iice], (alfa))

  return AVRG, Uprc, Lprc

# Create an array of day numbers with 0hr = init cond, 12 hr - daily means
# Assumed: runs start at 0 hr, if not - may need to change the logic 
# for finding the ic fields 
dnmbS = int(mtime.datenum([YR,MMS,DDS]))
dnmbE = int(mtime.datenum([YR,MME,DDE]))
RECS = [x + 0.5 for x in range(dnmbS, dnmbE + 1)]
RECS = np.array(RECS)

nexpts = len(ENMBS)
nrecs  = RECS.shape[0]
SMELT = np.zeros((nrecs,nexpts))
UPRC  = np.zeros((nrecs,nexpts))
LPRC  = np.zeros((nrecs,nexpts))

iens = -1
iprst = -1
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
      AI  = dcice['aice_d'].data.squeeze()
      AA    = dcice['melts_d'].data.squeeze()
      HS    = dcice['hs_d'].data.squeeze()

    AA = np.where(RMsk == 0, np.nan, AA)
    AI = np.where(RMsk == 0, np.nan, AI)
    HS = np.where(RMsk == 0, np.nan, HS)
    snow_mlt, uprc, lprc = spavrg_field(AA, HS, Acell)

    irec += 1
    print(f"Av. snow melt={snow_mlt:.3f} cm/day")

    SMELT[irec,iens] = snow_mlt
    UPRC[irec,iens]  = uprc
    LPRC[irec,iens]  = lprc

  
# Line colors:
import mod_gfs_cice_anls as mgfscice
importlib.reload(mgfscice)
CLRS = mgfscice.sens_tests_colors()

print("Plotting ...")

XT = RECS - np.floor(RECS[0])
xticks = np.arange(np.floor(XT[0]),np.ceil(XT[-1]))
yticks = np.arange(0.,2.,0.1)
dash = 0.05
show_prctl = False

sttl = f"Sp. avrg. snow melt cm/day datmUFS expts, {regn}\n {YR}/{MMS:02d}/{DDS:02d}-{YR}/{MME:02d}/{DDE:02d}"

plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.4, 0.8, 0.5])

LNS = []
for iens in range(nexpts):
  enmb = ENMBS[iens]
  smlt0 = SMELT[:,iens]
  uprc0 = UPRC[:,iens]
  lprc0 = LPRC[:,iens]
  clr0  = CLRS[iens,:]
  #line_lbl  = f"expt{enmb:02d}"
  line_lbl = mgfscice.sens_tests_info(enmb)
  ln1, = ax1.plot(XT, smlt0, 'o-', linewidth=2, color=clr0, label=line_lbl)
  LNS.append(ln1)

  nrc = len(smlt0)
  if show_prctl:
    for ik in range(nrc):
      yy1 = lprc0[ik]
      yy2 = uprc0[ik]
      xx0 = ik + 0.5
      xx1 = xx0 - dash
      xx2 = xx0 + dash
      ax1.plot([xx0,xx0], [yy1,yy2], '-', linewidth=1, color=clr0)
      ax1.plot([xx1,xx2], [yy1,yy1], '-', linewidth=1, color=clr0)
      ax1.plot([xx1,xx2], [yy2,yy2], '-', linewidth=1, color=clr0)
      

yl1 = 0
yl2 = np.nanmax(SMELT) * 1.05  
ax1.set_yticks(yticks)
ax1.set_xticks(xticks)
ax1.set_ylim(yl1, yl2)
ax1.grid('on')
ax1.set_xlabel('Forecast days')
ax1.set_title(sttl)

ax3 = plt.axes([0.1, 0.15, 0.6, 0.2])
lgd = plt.legend(handles=LNS, loc='upper left')
ax3.axis('off')

btx = 'meltsnow_expts_datmUFS.py'
bottom_text(btx, pos=[0.1,0.1])

