"""
  RMSE of ice conc. btw datm UFS experiments and NSIDC NRT fields

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

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help="hemisphere: north or south", type=str, required=True)
parser.add_argument("--yr", help=f"year of model run, default={YR}", type=int)
parser.add_argument("--ms", help=f"start month of data to plot, default={MM}", type=int)
parser.add_argument("--me", help=f"end month of NSIDCS data to plot, default={MM}", type=int)
parser.add_argument("--ds", help="Start: day in the start month to plot, default={DD}", type=int)
parser.add_argument("--de", help="End: day in the end month to plot, default={DD}", type=int)
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
DDE   = args.de if args.de else 16
ENMBS = args.enmb if args.enmb else None
plt_init = True  # show RMSE for init state if init. state file exists and saved by CICE6
  
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

def read_NSIDC(YR,MM,DD,regn,pthnsidc,varnm):
  if regn == 'south':
    flnsidc = f"sic_pss25_{YR}{MM:02d}{DD:02d}_am2_v06r00.nc"  
  else:
    flnsidc = f"sic_psn25_{YR}{MM:02d}{DD:02d}_am2_v06r00.nc"  

  with xarray.open_dataset(os.path.join(pthnsidc,flnsidc)) as ds_nsidc:
    A = ds_nsidc[varnm].data.squeeze()

  return A

pthdata = pths_ufs[node_nm]["MOM6"]["pthdata"]
pthnsidc = os.path.join(pthdata,f"NRT_NOAA_NSIDC_seaconc/{YR}")
if regn == 'south':
  RMsk = np.where(HH>=0, 0, 1)
  RMsk = np.where(hlat > -60., 0, RMsk)

# Create an array of day numbers with 0hr = init cond, 12 hr - daily means
dnmbS = int(mtime.datenum([YR,MMS,DDS]))
dnmbE = int(mtime.datenum([YR,MME,DDE]))
RECS = [dnmbS] + [x + 0.5 for x in range(dnmbS, dnmbE + 1)]
RECS = np.array(RECS)

nexpts = len(ENMBS)
nrecs  = RECS.shape[0]
RMSE = np.zeros((nrecs,nexpts))
iens = -1
for enmb in ENMBS:
  iens += 1
  irec = -1
  pthout_cice = pths_ufs[node_nm]["MOM6"]["pthcice"].format(enmb=enmb)

  for nn in range(nrecs):
    dnmb = RECS[nn]
    YR,MM,DD,hr = mtime.datevec(dnmb, round_hrs=True)[:4] 

    if hr == 0:
      # Initial state
      nsec0 = 0 
      flcice = f"iceh_ic.{YR}-{MM:02d}-{DD:02d}-{nsec0:05d}.nc"
      dflcice = os.path.join(pthout_cice,flcice)

      if not os.path.isfile(dflcice):
        print(f"Initial state file is missing, proceed without it ...")
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

    # Interpolated NSIDC obs fields:
    fliceout = f'NSIDC_iconc_interp_mesh025_{jdm}x{idm}_{YR}{MM:02d}_{regn}.nc'
    dfliceout = os.path.join(pthnsidc,fliceout)
    print(f'Loading interpolated ice conc {dfliceout}')
    with xarray.open_dataset(dfliceout) as dsint:
      AI = dsint['ice_conc'].isel(time=DD-1).squeeze()

    AA = np.where(RMsk == 0, np.nan, AA)
    AI = np.where(RMsk == 0, np.nan, AI)
    sqerr = (AA-AI)**2
   # Jice, Iice = np.where(~np.isnan(AA) & ~np.isnan(AI))
    Jice, Iice = np.where(((AA > 1.e-3) | (AI > 1.e-3)) & ~np.isnan(AA) & ~np.isnan(AI))

    Nice = len(Jice)
    assert(Nice>0), f"No ice grid found at {YR}/{MM:02d}/{DD:02d}"
    rmse_mo = np.sqrt(1./float(Nice)*np.sum(sqerr[Jice,Iice]))

    print(f"RMSE={rmse_mo:.2f}")
    irec += 1
    RMSE[irec,iens] = rmse_mo


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
sttl = f"RMSE btw iconc NSIDC and datmUFS expts, {YR}/{MMS:02d}/{DDS:02d}-{YR}/{MME:02d}/{DDE:02d}"

plt.ion()
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.1, 0.4, 0.8, 0.5])

LNS = []
for iens in range(nexpts):
  enmb = ENMBS[iens]
  rmse0 = RMSE[:,iens]
  clr0  = CLRS[iens,:]
  line_lbl  = f"expt{enmb:02d}"
  ln1, = ax1.plot(XT,rmse0, 'o-', linewidth=2, color=clr0, label=line_lbl)
  LNS.append(ln1)

ax1.set_xticks(XT)
ax1.grid('on')
ax1.set_xlabel('Forecast days')
ax1.set_title(sttl)

ax3 = plt.axes([0.1, 0.15, 0.6, 0.2])
lgd = plt.legend(handles=LNS, loc='upper left')
ax3.axis('off')

btx = 'calc_rmse_iconc_NSIDC.py'
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



