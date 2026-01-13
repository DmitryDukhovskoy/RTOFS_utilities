"""
  Plot ice concentration maps for Arctic / S. Ocean

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

expt = 'gfs_fcast'  # gfs - current f/cast in yaml
init_date = 20250714
init_hr = 6    # nominal hr, actual: -6 hrs for IAU, and -3 FHROT (f/cast hr rotation)
regn = 'south'

parser = argparse.ArgumentParser()
parser.add_argument("--regn", help=f"hemisphere: north or south, default={regn}", type=str)
parser.add_argument("--init", help=f"init date, default={init_date}", type=int)
parser.add_argument("--ihr", help=f"init hour, default={init_hr}", type=int)
parser.add_argument("--fhr", help=f"forecast hour to plot: 6, 12, ...,390, =0 - init. cond.", type=int, required=True)
parser.add_argument("--enmb", help="experiment number: 1, 2, ...", type=int)
args = parser.parse_args()
  
enmb      = args.enmb if args.enmb else None
init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
fhr       = args.fhr if args.fhr is not None else None
regn      = args.regn if args.regn else regn

dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
YRI,MMI,DDI,hrI = mtime.datevec(dnmbI, round_hrs=True)[:4]

dnmb0 = dnmbI + fhr/24
yr0,mm0,dd0,hr0 = mtime.datevec(dnmb0, round_hrs=True)[:4]
YR,MM,DD = mtime.datevec(dnmb0)[:3]
  
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
if expt == 'gfs_fcast':
  pthgrid = pths_ufs[node_nm]["GFS"]["pthgrid"]
  pthcomroot = pths_ufs[node_nm]["GFS"]["pthcomroot"]
  subdir_ice = pths_ufs[node_nm]["GFS"]["subdir_ice"].format(YR=YRI,MM=MMI,DD=DDI,hr=hrI)
  pthoutp = os.path.join(pthcomroot,subdir_ice)
else:
  pthgrid = pths_ufs[node_nm]["MOM6"]["pthgrid"]
  pthoutp = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/ufs_datm_mx025/expt{enmb:02d}/cice6"

dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")
    
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

jdm, idm = HH.shape

# Get date:
plot_init = fhr == 0  # initial conditions


if expt == 'gfs_fcast':
  if plot_init:
    flinp = "gfs.t06z.ic.nc"
  else:
    flinp = f"gfs.t06z.6hr_avg.f{fhr:03d}.nc"
  varnm = 'Tsfc_h'
else:
  if plot_init:
    flinp = f"iceh_ic.{yr0}-{mm0:02d}-{dd0:02d}-{nsec0:05d}.nc"
  else:
    flinp = f"iceh.{yr0}-{mm0:02d}-{dd0:02d}.nc"
  varnm = 'Tsfc_d'

dflice = os.path.join(pthoutp,flinp)

print(f"Processing {YR}/{MM}/{DD} {hr0:02d}:00, {expt} init {YRI}/{MMI:02d}/{DDI:02d} {hrI:02d}:00\n{dflice}")
with xarray.open_dataset(dflice) as dcice:
  AA = dcice[varnm].data.squeeze()
  IConc = dcice['aice_h'].data.squeeze()

# Mask land and no ice regions
AA = np.where(HH >= 0, np.nan, AA)
AA[IConc <= 1.e-8] = 100.

plt.ion()


clrmp = mclrmps.colormap_cold_warm(ins_white=False)
rmin = -30.
rmax = 10.
clrmp.set_bad(color=[0.2, 0.2, 0.2])
clrmp.set_over(color=[1,1,1])

if regn == 'south':
  m = Basemap(projection='spstere',boundinglat=-50,lon_0=180,resolution='l')
#lons, lats = m.makegrid(idim, jdim) # get lat/lons of ny by nx evenly spaced grid.
parallels = np.arange(-80,-10,10.)
meridians = np.arange(-360,359.,45.)
xl1 = -8.6e6
xl2 = -0.8e6
yl1 = xl1
yl2 = xl2 

xh, yh = m(hlon,hlat) # GFS coords
fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.08, 0.13, 0.83, 0.83])
m.drawparallels(parallels,labels=[0,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10)
img1 = ax1.pcolormesh(xh,yh,AA, cmap=clrmp, vmin=rmin, vmax=rmax)
#img1 = ax1.pcolormesh(xh,yh,sqerr, cmap=clrmp, vmin=rmin, vmax=rmax)
#img1 = ax1.pcolormesh(xh,yh,np.abs(AA-AI), cmap=clrmp, vmin=rmin, vmax=rmax)
Tmin_ant = np.nanmin(AA[hlat<-45])
Tmin_arc = np.nanmin(AA[hlat>55])
if regn == 'south':
  Tmin = Tmin_ant
else:
  Tmin = Tmin_arc

ax1.set_title(f'{expt} init {init_date}/{init_hr:02d}, {varnm}, fcast {fhr}hr  {YR}/{MM:02d}/{DD:02d} {hr0:02d}, Tmin={Tmin:0.2f}')
ax1.set_xlim([xl1, xl2]) 
ax1.set_ylim([yl1, yl2]) 
ax1.invert_yaxis()
ax1.invert_xaxis()


# Colorbars
ax2 = fig1.add_axes([0.15, 0.1, 0.7, 0.02])
clb = plt.colorbar(img1, cax=ax2, orientation='horizontal', extend='both')
ax2.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_xticklabels(ax2.get_xticks())
ticklabs = clb.ax.get_xticklabels()
clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)


sinfo = pthoutp + '\n'
sinfo = sinfo + flinp
ax3 = fig1.add_axes([0.02, 0.03, 0.8, 0.06])
ax3.text(0, 0, sinfo, fontsize=8)
ax3.axis('off')


btx = 'maps_surfaceT_GFS.py'
bottom_text(btx, pos=[0.1,0.01], fsz=8)



