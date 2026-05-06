"""
  RMSE of ice conc. btw datm UFS experiments, RTOFS CICE4  and NSIDC NRT fields

  NSIDC fields from 
  https://noaadata.apps.nsidc.org/NOAA/G02202_V6/north/daily/2025/

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib  
import xarray
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

init_hr = 0

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere: north or south", type=str, required=True)
parser.add_argument("--init", help="Init date of the f/cast YYYYMMDD", 
                    choices=[20250704, 20241231], required=True, type=int)
parser.add_argument("--ndays", help="N of fcast days to calc. RMSE", required=True, type=int)
parser.add_argument(
    "--enmb",
    help="List of DATM UFS experiments that start on init",
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
init_date = args.init
ndays     = args.ndays
ENMBS     = args.enmb if args.enmb else None
PRST      = args.prst if args.prst else []
plt_init = True  # show RMSE for init state if init. state file exists and saved by CICE6

fhr = ndays // 24

# Error in NSIDC ice concentration fields Northern h/sphere:
if regn == 'north':
  NSIDC_err = mtime.datenum_v2([[2024,7,12],[2025,7,27]], ref_day0=False)
else:
  NSIDC_err = None
  
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
#pthnsidc = os.path.join(pthdata,f"NRT_NOAA_NSIDC_seaconc/{YR}")

RMsk = np.where(HH>=0, 0, 1)
if regn == 'south':
  RMsk = np.where(hlat > -60., 0, RMsk)
elif regn == 'north':
  RMsk = np.where(hlat < 50, 0, RMsk)

def rmse2d(AA,AI):
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

def rmse_NSIDC(AA, pthdata, regn, dnmb):
  YR, MM, DD  = mtime.datevec(dnmb)[:3]
  if int(dnmb) in NSIDC_err: 
    # Error ice conc fields
    rmse_mo = np.nan
  else:
    # Interpolated NSIDC obs fields:
    pthnsidc = os.path.join(pthdata,f"NRT_NOAA_NSIDC_seaconc/{YR}")
    fliceout = f'NSIDC_iconc_interp_mesh025_1080x1440_{YR}{MM:02d}_{regn}.nc'
    dfliceout = os.path.join(pthnsidc,fliceout)
    print(f'Loading interpolated ice conc {dfliceout}')
    with xarray.open_dataset(dfliceout) as dsint:
      AI = dsint['ice_conc'].isel(time=DD-1).squeeze()

    AA = np.where(RMsk == 0, np.nan, AA)
    AI = np.where(RMsk == 0, np.nan, AI)
    rmse_mo = rmse2d(AA,AI)

  return rmse_mo

# Create an array of day numbers with 0hr = init cond, 12 hr - daily means
# Assumed: runs start at 0 hr, if not - may need to change the logic 
# for finding the ic fields 
# Note: RTOFS cice output are instanteneous at the end of 24-hr cycles
dnmbS = int(np.floor(mtime.rdate2datenum(init_date*100+init_hr)))
YR, MMS, DDS = mtime.datevec(dnmbS)[:3]
dnmbE = int(dnmbS + ndays)
YRE, MME, DDE = mtime.datevec(dnmbE)[:3]
RECS = [dnmbS] + [x + 0.5 for x in range(dnmbS, dnmbE + 1)]
RECS = np.array(RECS)
RECS_rtofs = np.arange(dnmbS,dnmbS+np.min([ndays,8])+1)  # time at 24:00, 8 days forecasts only

nprst  = len(PRST)
nexpts = len(ENMBS)
nrecs  = RECS.shape[0]
RMSE = np.zeros((nrecs,nexpts))
if nprst > 0:
  RMSEp = np.zeros((nrecs, nprst))
Aprst = None

iens = -1
iprst = -1
# Run datm UFS first, then RTOFS
for enmb in ENMBS:
  iens += 1
  irec = -1
  pthout_cice = pths_ufs[node_nm]["MOM6"]["pthcice"].format(enmb=enmb)

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

    rmse_mo = rmse_NSIDC(AA, pthdata, regn, dnmb)

    irec += 1
    if track_prst:
      if Aprst is  None:
        Aprst = AA.copy()
        iprst += 1
      rmse_prst = rmse_NSIDC(Aprst, pthdata, regn, dnmb)
      print(f"RMSE={rmse_mo:.3f}  RMSE_prst={rmse_prst:.3f}")
      RMSEp[irec,iprst] = rmse_prst
    else:
      print(f"RMSE={rmse_mo:.3f}")

    RMSE[irec,iens] = rmse_mo

  Aprst = None
 
# RTOFS f/cast:
if dnmbS < mtime.datenum([2025,8,1]):
  rtofs_vers = "2.4"
else:
  rtofs_vers = "2.5"
 
pth_rtofs = f'/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/RTOFS/cice4_fcast/{init_date}/interp_mesh025'
irec = -1
RMSErtofs = np.zeros((len(RECS_rtofs)))
for dnmb in RECS_rtofs:
  YR,MM,DD,hr = mtime.datevec(dnmb, round_hrs=True)[:4]  
  fhr0 = int((dnmb - dnmbS) * 24)
  flname = f'iconc_RTOFSv{rtofs_vers}_{init_date}_fhr{fhr0:03d}_1080x1440.nc'
  dflname = os.path.join(pth_rtofs, flname)

  print(f"Reading RTOFS iconc: {dflname}")
  with xarray.open_dataset(dflname) as dsrtofs:
    AA = dsrtofs['ice_conc'].isel(time=0).squeeze()

  rmse_mo = rmse_NSIDC(AA, pthdata, regn, dnmb)
  irec += 1
  RMSErtofs[irec] = rmse_mo 


# Line colors:
import mod_gfs_cice_anls as mgfscice
importlib.reload(mgfscice)
CLRS = mgfscice.sens_tests_colors()

clr_rtofs = [0.,0.,0.]

print("Plotting ...")

XT = RECS - np.floor(RECS[0])
XTr = RECS_rtofs - np.floor(RECS[0])  
xticks = np.arange(np.floor(XT[0]),np.ceil(XT[-1]+1))
yticks = np.arange(0.,1.,0.05)
sttl = f"RMSE iconc btw datmUFS, RTOFSv{rtofs_vers} vs NSIDC, {regn}\n"
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
  line_lbl = mgfscice.sens_tests_info(enmb)
  ln1, = ax1.plot(XT,rmse0, 'o-', linewidth=2, color=clr0, label=line_lbl)
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

# RTOFS f/cast:
line_lbl = "RTOFS"
ln1, = ax1.plot(XTr, RMSErtofs, 'o-', linewidth=2, color=clr_rtofs, label=line_lbl)
LNS.append(ln1)

yl1 = 0
yl2 = np.nanmax(RMSE) * 1.05 
#yl2 = 0.5

if nprst > 0:
  ylP = np.nanmax(RMSEp) * 1.05
  yl2 = np.max([yl2, ylP])
 
ax1.set_yticks(yticks)
ax1.set_xticks(xticks)
ax1.set_ylim(yl1, yl2)
ax1.set_xlim(0,xticks[-1])
ax1.grid('on')
ax1.set_xlabel('Forecast days')
ax1.set_title(sttl)

ax3 = plt.axes([0.1, 0.15, 0.6, 0.2])
lgd = plt.legend(handles=LNS, loc='upper left')
ax3.axis('off')

btx = 'calc_rmse_iconc_datmUFS_RTOFS.py'
bottom_text(btx, pos=[0.1,0.1])

f_chck = False
if f_chck:
  clrmp = mclrmps.colormap_conc()
  rmin = 0.
  rmax = 1.
  clrmp.set_bad(color=[0.2, 0.2, 0.2])

  clrmp_dlt = mclrmps.colormap_uv()
  dmin = -1
  dmax = 1
  clrmp_dlt.set_bad(color=[0.2, 0.2, 0.2])

  if regn == 'south':
    m = Basemap(projection='spstere',boundinglat=-50,lon_0=180,resolution='l')
  #lons, lats = m.makegrid(idim, jdim) # get lat/lons of ny by nx evenly spaced grid.
  parallels = np.arange(-80,-10,10.)
  meridians = np.arange(-360,359.,45.)
  xl1 = -8.e6
  xl2 = -1.2e6
  yl1 = xl1
  yl2 = xl2 
 
  xh, yh = m(hlon,hlat) # GFS coords
  plt.clf()
  ax1 = plt.axes([0.05, 0.55, 0.4, 0.4])
  m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
  m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
  img1 = ax1.pcolormesh(xh,yh,AA, cmap=clrmp, vmin=rmin, vmax=rmax)
  #img1 = ax1.pcolormesh(xh,yh,sqerr, cmap=clrmp, vmin=rmin, vmax=rmax)
  #img1 = ax1.pcolormesh(xh,yh,np.abs(AA-AI), cmap=clrmp, vmin=rmin, vmax=rmax)
  ax1.set_title(f'UFS CICE iconc {YR}/{MM:02d}/{DD:02d}')
  ax1.set_xlim([xl1, xl2]) 
  ax1.set_ylim([yl1, yl2]) 
  ax1.invert_yaxis()
  ax1.invert_xaxis()

  # Interpolated iconc
  ax2 = plt.axes([0.55, 0.55, 0.4, 0.4])
  m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
  m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
  img2 = ax2.pcolormesh(xh, yh, AI, cmap=clrmp, vmin=rmin, vmax=rmax)
  ax2.set_title('NSIDC iconc interp to mesh025')
  ax2.set_xlim([xl1, xl2])      
  ax2.set_ylim([yl1, yl2])      
  ax2.invert_yaxis()
  ax2.invert_xaxis()

  ax21 = plt.axes([0.05, 0.1, 0.4, 0.4])
  m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
  m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
  ax21.pcolormesh(xh,yh,np.abs(AA-AI), cmap=clrmp, vmin=rmin, vmax=rmax)
  ax21.set_title(f'|err| UFS CICE vs NSIDC {YR}/{MM:02d}/{DD:02d}')
  ax21.set_xlim([xl1, xl2])      
  ax21.set_ylim([yl1, yl2])      
  ax21.invert_yaxis()
  ax21.invert_xaxis()
  
  ax22 = plt.axes([0.55, 0.1, 0.4, 0.4])
  m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
  m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
  img2 = ax22.pcolormesh(xh,yh,(AA-AI), cmap=clrmp_dlt, vmin=dmin, vmax=dmax)
  ax22.set_title(f'diff UFS CICE vs NSIDC {YR}/{MM:02d}/{DD:02d}')
  ax22.set_xlim([xl1, xl2])      
  ax22.set_ylim([yl1, yl2])      
  ax22.invert_yaxis()
  ax22.invert_xaxis()

  
  # Colorbars
  ax3 = fig1.add_axes([0.05, 0.05, 0.4, 0.02])
  clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
  ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax3.set_xticklabels(ax3.get_xticks())
  ticklabs = clb.ax.get_xticklabels()
  clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  ax4 = fig1.add_axes([0.55, 0.05, 0.4, 0.02])
  clb = plt.colorbar(img2, cax=ax4, orientation='horizontal', extend='max')
  ax4.xaxis.set_ticks(list(np.linspace(dmin,dmax,11)))
  ax4.set_xticklabels(ax4.get_xticks())
  ticklabs = clb.ax.get_xticklabels()
  clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  bottom_text(btx, pos=[0.1,0.01])



