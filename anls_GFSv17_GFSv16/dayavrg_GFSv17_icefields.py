"""

  Produce daily averaged from 6hr f/casts

  in CICE6 output:
  albedo is in percent
  averaged:  "averaged for coszen>0, weighted by aice" ;
  so for polar night albedo = 0
 
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
import pandas as pd
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
import mod_regmom as mrmom 
import mod_sis2_relax as msisrlx

# fcast init date:
init_date = 20251231
init_hr = 0
regn = 'global'
fhrS = 6    # forecast start hour for averaging
fhrE = 384  # forecast end hour
dhr  = 6    # output freq., hours
hrav = 24   # averaging period, hrs

parser = argparse.ArgumentParser()
parser.add_argument("--field", help=f"Field to interpolate: snow thkciness or ice thickn",
                    choices=['hsnow', 'ithkn', 'snwed', 'iconc', 'albd'], required=True, type=str)
parser.add_argument("--init", help=f"init date, default {init_date}", type=int)
parser.add_argument("--init_hr", help=f"init hour, default {init_hr}", type=int)
parser.add_argument("--fhrS", help=f"start: forecasts hour, default={fhrS}", type=int)
parser.add_argument("--fhrE", help=f"end: forecasts hour, default={fhrE}", type=int)
parser.add_argument("--dhr", help=f"Output freq., hours:  default={dhr}", type=int)
parser.add_argument("--hrav", help=f"Averaging period, hours", type=int)
args = parser.parse_args()
  
field_name = args.field if args.field else None
init_date  = args.init  if args.init else init_date
init_hr    = args.init_hr   if args.init_hr is not None else init_hr
fhrS       = args.fhrS  if args.fhrS is not None else fhrS
fhrE       = args.fhrE  if args.fhrE is not None else fhrE
dhr        = args.dhr   if args.dhr is not None else dhr
hrav       = args.hrav  if args.hrav is not None else hrav

# Init date:
dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
YR, MM, DD = mtime.datevec(dnmbI)[:3]


syst_info = os.uname() 
machine = syst_info.nodename
  
if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea compute node:", machine)
  node_nm = "gaea"
elif 'an' in machine:
  print("Running on PPAN node:", machine)
  node_nm = "ppan"
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
# Mask of polare regions wihtout land 
latS = -50.
latN = 50.
LMsk = HH < 0
LMsk[(hlat > latS) & (hlat < latN)] = False


match field_name:
  case 'hsnow':
    varnm = 'hs_h'
    attr_str_long = 'grid cell mean snow thickness'
    attr_units = 'm'
  case 'ithkn':
    varnm = 'hi_h'
    attr_str_long = 'grid cell mean ice thickness'
    attr_units = 'm'
  case 'iconc':
    varnm = 'aice_h'
    attr_str_long = 'sea ice partial area'
    attr_units = 'partial coverage'
  case 'albd':
    varnm = 'albsni_h'   # prcnt
    attr_str_long = 'snow/ice surface albedo'
    attr_units = 'percent'
  case _:
    raise VarnameError(f"{field_name} not recognized")


# Time index
TMF = np.array([x for x in range(fhrS, fhrE+dhr, dhr)])
nrec_avrg = hrav // dhr

# Find N of records in snow/ice file
idate_str = f"{YR}{MM:02d}{DD:02d}"
init_hr_str = f"{init_hr:02d}"
pthgfs = '/gpfs/f6/gfs-cpu/scratch/Xiao.Luo/DIAG/RETROV17/'
pthfcast = os.path.join(pthgfs,f"gfs.{idate_str}",init_hr_str,'model','ice','history')


# Averaged fields:
pthavrg = f"/gpfs/f6/sfs-emc/proj-shared/Dmitry.Dukhovskoy/GFSv17/{idate_str}/daily"
os.makedirs(pthavrg, exist_ok=True)

assert hrav == 24, f"Daily averaging assumed, hrav={hrav} hrs"

irec_avrg = 0
iday = 0
A3d = None
A3d_list = []
TM = []
for FH in TMF:
  flgfs = f"gfs.t00z.6hr_avg.f{FH:03d}.nc"
  dflgfs = os.path.join(pthfcast, flgfs)
  print(f"Reading GFSv17: {dflgfs}")

  with xarray.open_dataset(dflgfs) as ds_gfs:
    A2d = ds_gfs[varnm].values.squeeze()

  if irec_avrg == 0:
    Asum = A2d.copy() 
  else:
    Asum += A2d

  irec_avrg += 1
  if irec_avrg == nrec_avrg:
    # assumed 24 hr averaging !
    A2d = Asum / irec_avrg
    iday += 1
    dnmb = mtime.datenum([YR,MM,DD]) + iday 

    A3d_list.append(A2d)
    TM.append(dnmb)

    yr0, mm0, dd0 = mtime.datevec(dnmb)[:3]
    print(f"Averging {field_name} {yr0}/{mm0:02d}/{dd0:02d}")

    irec_avrg = 0


# Finished, save netcdf
A3d = np.stack(A3d_list, axis=0)
dnmb_ref = mtime.datenum([1900, 1, 1])
DV = mtime.datevec(dnmb_ref)[:3]
TM = (np.array(TM) - dnmb_ref - 0.5).astype("float64") # daily mean

assert A3d.shape[1:] == hlon.shape, "Check A3d shape [time, jdim, idim]?"

jdim, idim = hlon.shape
darr_cice = xarray.DataArray(A3d, dims=("time","jdim","idim"),\
                     coords={"time": TM,\
                             "jdim": np.arange(jdim),\
                             "idim": np.arange(idim)})
dset = xarray.Dataset({varnm: darr_cice})
dset[varnm].attrs['long_name'] = attr_str_long
dset[varnm].attrs['units'] = attr_units
dset["time"].attrs = {
     "long_name": "time",
     "units": f"days since {DV[0]}-{DV[1]:02d}-{DV[2]:02d}",
}

# Add global attributes:
dset.attrs['title']       = f'GFSv17 forecast {attr_str_long} daily averaged'
dset.attrs['source']      = 'dayavrg_GFSv17_icefields.py'

ndays = int(np.ceil((fhrE-fhrS)/24))
dayS = fhrS // 24
dayE = fhrE // 24
flout = f"gfsv17_{field_name}_init{YR}{MM:02d}{DD:02d}_days{dayS:03d}_{dayE:03d}.nc"
dfliceout = os.path.join(pthavrg, flout)
print(f'Dumping averaged {field_name} --> {dfliceout}\n')
dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')


f_chck = False
if f_chck:
  regn = 'south' 
  plt_intrp = True  
  if field_name == 'iconc':
    clrmp = mclrmps.colormap_conc()
    rmin = 0.
    rmax = 1. 
  elif field_name == 'ithkn':
    clrmp = mclrmps.colormap_ice_thkn()
    rmin = 0.
    rmax = 3.

  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  clrmp.set_under(color=[1,1,1])

  if regn == 'north':
    m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
    parallels = np.arange(40,89,10.)
    meridians = np.arange(-360,359.,45.)
    # Subset to avoid inf for opposite pole in stereogr proj
    jcut = np.min(np.where(hlat >= 0)[0])
    LAT0 = hlat[jcut:,:]
    LON0 = hlon[jcut:,:]
    A0 = HSint[jcut:,:]
  elif regn == 'south':
    m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
    parallels = np.arange(-80,0,10.)
    meridians = np.arange(-360,359.,45.)
    # Subset to avoid inf for opposite pole in stereogr proj
    jcut = np.min(np.where(LAT >= 0)[0])
    LAT0 = hlat[:jcut,:]
    LON0 = hlon[:jcut,:]
    A0 = HSint[:jcut,:]

  xh, yh = m(LON0, LAT0)              # mesh025 GFSv17


  plt.ion() 
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.12, 0.12, 0.8, 0.8])

  ax1.cla()
  m.drawcoastlines()
  m.drawparallels(parallels,labels=[0,0,0,0])
  m.drawmeridians(meridians,labels=[0,0,0,0])

  img = ax1.pcolormesh(xh, yh, A0, cmap=clrmp, vmin=rmin, vmax=rmax)
  sttl = f"{field_name} GFSv17 interp mesh025 \n{flgfs}"

  ax1.set_title(sttl)

  ax2 = fig1.add_axes([0.9,0.1,0.015,0.8])
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
  ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=10)
        
  btx = 'dayavrg_GFSv17_icefields.py'
  bottom_text(btx) 



