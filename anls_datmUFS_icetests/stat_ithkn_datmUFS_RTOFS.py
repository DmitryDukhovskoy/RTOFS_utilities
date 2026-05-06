"""
  Compare statistics of ice thickness in datmUFS with RTOFS CICE5 IC vs 
  RTOFSv2.4/2.5 forecasts

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
parser.add_argument("--regn", help="hemisphere: north or south", 
                    choices=['north','south'], type=str, required=True)
parser.add_argument("--init", help="Init date of the f/cast YYYYMMDD", 
                    choices=[20250704, 20241231], required=True, type=int)
parser.add_argument("--fday", help="Fcast day: 0,..,8", required=True, type=int)
parser.add_argument(
    "--enmb",
    help="DATM UFS experiment that start on init",
    type=int,
    required=True
)
args = parser.parse_args()
  
regn      = args.regn if args.regn else None
init_date = args.init
fday      = args.fday
enmb      = args.enmb if args.enmb else None
plt_init = True  # show RMSE for init state if init. state file exists and saved by CICE6

fhr = fday // 24

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


# Note: RTOFS cice output are instanteneous at the end of 24-hr cycles
# Get date: 
plot_init = fday == 0  # initial conditions for datmUFS

dnmbI = int(np.floor(mtime.rdate2datenum(init_date*100+init_hr)))
YRI, MMI, DDI = mtime.datevec(dnmbI)[:3]

dnmbF = int(dnmbI + fday)
YRF, MMF, DDF = mtime.datevec(dnmbF)[:3]


iens = -1
iprst = -1
# Run datm UFS first, then RTOFS
pthoutp = pths_ufs[node_nm]["MOM6"]["pthcice"].format(enmb=enmb)

if plot_init:
  nsec0 = 0
  flinp = f"iceh_ic.{YRI}-{MMI:02d}-{DDI:02d}-{nsec0:05d}.nc"
else:
  flinp = f"iceh.{YRF}-{MMF:02d}-{DDF:02d}.nc"
varnm = 'hi_d'

dflice = os.path.join(pthoutp, flinp)
print(f"Processing {YRF}/{MMF}/{DDF}, expt {enmb:02d} init {YRI}/{MMI:02d}/{DDI:02d}\n{dflice}")
with xarray.open_dataset(dflice) as dcice:
  HIufs = dcice['hi_d'].data.squeeze()
  AIufs = dcice['aice_d'].data.squeeze()

# RTOFS f/cast:
if dnmbI < mtime.datenum([2025,8,1]):
  rtofs_vers = "2.4"
else:
  rtofs_vers = "2.5"
 
pth_rtofs = f'/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/RTOFS/cice4_fcast/{init_date}/interp_mesh025'
flname = f'ithkn_RTOFSv{rtofs_vers}_{init_date}_fhr{fhr:03d}_1080x1440.nc'
dflname = os.path.join(pth_rtofs, flname)
print(f"Reading RTOFS iconc: {dflname}")
with xarray.open_dataset(dflname) as dsrtofs:
  HIrtofs = dsrtofs['ice_thkn'].isel(time=0).data.squeeze()

flname = f'iconc_RTOFSv{rtofs_vers}_{init_date}_fhr{fhr:03d}_1080x1440.nc'
dfrtofs = os.path.join(pth_rtofs, flname)
print(f"Reading RTOFS iconc: {dfrtofs}")
with xarray.open_dataset(dfrtofs) as dsrtofs:
  AIrtofs = dsrtofs['ice_conc'].isel(time=0).data.squeeze()

# Applying masks to isolate study region:
if regn == 'north':
  RGmsk = (HH < 0) & (hlat > 50.)
elif regn == 'south':
  RGmsk = (HH < 0) & (hlat < -50.)

AIufs[~RGmsk] = np.nan
HIufs[~RGmsk] = np.nan

AIrtofs[~RGmsk] = np.nan
HIrtofs[~RGmsk] = np.nan

# Discard open water 
# Consider only grid points with ice in both cases
puny = 1.e-2   # exclude small concentrations cases
IMsk_ufs   = (AIufs > puny) & (np.isfinite(AIufs))
IMsk_rtofs = (AIrtofs > puny) & (np.isfinite(AIrtofs))

HIufs[~IMsk_ufs | ~IMsk_rtofs] = np.nan
HIrtofs[~IMsk_ufs | ~IMsk_rtofs] = np.nan

Msk_valid = np.isfinite(HIufs) & np.isfinite(HIrtofs)
Nvalid = np.sum(Msk_valid)

# Pre-extract valid grid points:
HIu = HIufs[Msk_valid]
HIr = HIrtofs[Msk_valid]

# Calc abs err:
abs_err = np.abs(HIu - HIr) 
diff_err = (HIu - HIr)
rmse = np.sqrt(1/Nvalid * np.sum( (HIu - HIr)**2 ))

hbin = np.append(np.arange(0, 4.5, 0.5), 1000.)


# Histogram counts:
hist_ufs, _   = np.histogram(HIu, bins=hbin)
hist_rtofs, _ = np.histogram(HIr, bins=hbin)

# Convert to frequency:
hist_ufs = hist_ufs / np.sum(hist_ufs)
hist_rtofs = hist_rtofs / np.sum(hist_rtofs)


# Bin-average statistics:
Ibin = np.digitize(HIu, hbin)
nbins = len(hbin) - 1
err_bin = np.full(nbins, np.nan)
rmse_bin = np.full(nbins, np.nan)

for ii in range(1, nbins+1):
  msk = Ibin == ii
  nbin = np.sum(msk)
  if np.any(msk):
    err_bin[ii-1] = np.mean(diff_err[msk])
    rmse_bin[ii-1] = np.sqrt(np.mean( (HIu[msk] - HIr[msk])**2 )) 


# Line colors:
import mod_gfs_cice_anls as mgfscice
importlib.reload(mgfscice)
CLRS = mgfscice.sens_tests_colors()

clr_rtofs = [0.,0.,0.]

print("Plotting ...")

line_lbl = mgfscice.sens_tests_info(enmb)
sinfo = f"datmUFS expt{enmb:02d} {line_lbl}, RTOFSv{rtofs_vers} INIT: {init_date} {regn}\n"
sinfo = sinfo + f"datmUFS: {dflice}\n"
sinfo = sinfo + f"RTOFS: {dfrtofs}"

plt.ion()
fig1 = plt.figure(1,figsize=(9,9))

# Hist:
bin_cntr = 0.5 * (hbin[:-1] + hbin[1:])
bin_dw  = np.diff(hbin)
bin_dw[-1] = bin_dw[-2]
bin_cntr[-1] = hbin[-2] + bin_dw[-1]/2 
hist_xticks = hbin
hist_xticks[-1] = hist_xticks[-2] + bin_dw[-1]

# shrink a bit so bars do not overlap
bw = 0.4 * bin_dw  

sttl1 = f"ithkn, {init_date}, FDAY={YRF}/{MMF}/{DDF}, {regn}"
clr_ufs = [0.2, 0.5, 0.8]
clr_rtofs = [0.8, 0.3, 0.3]
clr_rmse = [1,0.6,0.]
clr_err = [0.,0.4,0.9]

plt.clf()
ax1 = plt.axes([0.08, 0.6, 0.42, 0.35])
ax1.bar(bin_cntr - bw/2, hist_ufs,
        width=bw, color=clr_ufs, label='DATM UFS')

ax1.bar(bin_cntr + bw/2, hist_rtofs,
        width=bw, color=clr_rtofs, label=f'RTOFSv{rtofs_vers}')

ax1.set_xlabel('Ice thickness (m)')
ax1.set_ylabel('Frequency')
ax1.xaxis.set_ticks(hist_xticks)
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_title(sttl1)

sttl2 = f"RMSE ithkn RTOFS vs datmUFS, {regn}"
ax2 = plt.axes([0.57, 0.6, 0.42, 0.35])
ax2.bar(bin_cntr, rmse_bin, width=0.45, color=clr_rmse)
ax2.set_xlabel('Ice thickness (m)')
ax2.xaxis.set_ticks(hist_xticks)
ax2.grid(True, alpha=0.3)
ax2.set_title(sttl2)

sttl3 = f"Error ithkn datmUFS-RTOFS, {regn}"
ax3 = plt.axes([0.08, 0.15, 0.42, 0.35])
ax3.bar(bin_cntr, err_bin, width=0.45, color=clr_err)
ax3.set_xlabel('Ice thickness (m)')
ax3.xaxis.set_ticks(hist_xticks)
ax3.grid(True, alpha=0.3)
ax3.set_title(sttl3)

ax5 = plt.axes([0.02, 0.04, 0.8, 0.05])
ax5.text(0,0, sinfo, fontsize=8)
ax5.axis('off')

btx = 'stat_ithkn_datmUFS_RTOFS.py'
bottom_text(btx, pos=[0.05,0.02])



