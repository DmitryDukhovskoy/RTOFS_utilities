"""
  General prediction using ML model
  and any input fields

  Input fields (predictors) are from GFSv17 analysis
  GDAS atm
  JEDI SOCA ice / ocean

  Specify ice concentration fields used as an input for ML emulator

  Random Forest or 
  Hist Gradient Boost Regressor (decision trees) predictor

  Run N days predictions
  predictions will be performed from sdate to edate
  using available ERA5 daily SAT fields

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import xarray as xr
import matplotlib.colors as colors
import argparse
from yaml import safe_load
import joblib
import json

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
    os.path.join(PPTHN, 'MyPython'),
    os.path.join(PPTHN, 'MyPython', 'mom6_utils')
])
from mod_utils_fig import bottom_text, colorbar_horiz
import mod_time as mtime
import mod_mom6 as mmom6
import mod_ml_emulator as mml

# See mod_icepredict.py for more info about the models
#    Random Forest (all 1993-2025):
#      RF1:  Ntrees 100
#      RF2:  Ntrees 200
#    Hist. Grad Boost Regressor (decision tree)
#      GBR1: max_leaf = 8,  max_iter=500,  max_depth=8
#      GBR2: max_leaf = 31, max_iter=1000, max_depth=10
#      GBR3: max_leaf = 63, max_iter=1500, max_depth=15
parser = argparse.ArgumentParser()
parser.add_argument("--model", help="ML model to use",
                   choices=['rf1','rf2','gbr2','gbr3'],
                   type=str,
                   required=True)
parser.add_argument("--sdate", help="Start prediction date YYYYMMDD", required=True, type=int)
parser.add_argument("--edate", help="End prediction date YYYYMMDD, skip if edate=sdate",
                    type=int)
parser.add_argument("--regn", help="Region to process", choices=['north','south'],
                    required=True, type=str)
parser.add_argument("--sinput", help="=1: Save created input fields for inspection (default)",
                    choices=[0,1], type=int)
parser.add_argument("--save", 
  help="Save predicted ithkn: 0=no (default), 1=npz (grid pnts), 2=netcdf (whole domain))",
  choices=[0,1,2], 
  default=0, 
  type=int)
args = parser.parse_args()

ml_model  = args.model
sdate     = args.sdate
edate     = args.edate if args.edate is not None else args.sdate
regn      = args.regn
save_input = args.sinput == 1    # save created input (predictor) fields
save_npz   = args.save == 1      # save numpy binary
save_nc    = args.save == 2      # save netcdf

if save_npz:
  print("Output will be saved to numpy binary npz")
elif save_nc:
  print("Output will be saved to netcdf")
else:
  print(" \n!!!!   PREDICTIONS WILL NOT BE SAVED !!!\n")

f_plt = False
standz = False    # standardized predictors only for lin. regr
iconc_min = 0.05  # Discard too low ice conc. predictions
sst_max = 5.      # Discard ice in too warm ocean

fyaml = "paths_ML.yaml"
with open(fyaml) as ff:
  config_ml = safe_load(ff)

# Requested start / end time - those may change
# based on saved ERA5 time
dnmbS = mtime.rdate2datenum(sdate)
YRS, MMS, DDS = mtime.datevec(dnmbS)[:3]
dnmbE = mtime.rdate2datenum(edate)
YRE, MME, DDE = mtime.datevec(dnmbE)[:3]

# ML emulator lat boundaries:
lat0 = config_ml[regn]["lat0"]

# Load model parameters and RF model:
model_name = config_ml["ml"][ml_model]["mlname"].format(regn=regn)
pthmodel = config_ml["ml"][ml_model]["pthmodel"]
model_file = os.path.join(pthmodel, model_name + ".pkl")
print(f"Reading {ml_model} object from {model_file}")
mlem = joblib.load(model_file)

# Training period:
info_file = os.path.join(pthmodel, model_name + "_info.json")
with open(info_file, "r") as f:
  info = json.load(f)

YS           = info["training_yrS"]
YE           = info["training_yrE"]
dxy          = info["dxy"]             # ice length scale
Tfrz         = info["Tfrz"]
sqrt_frzdays = info["sqrt_frzdays"]
intgr_time   = info["intgr_time"]
ndays_era    = info["ndays_era"]


# ML Predictors / Input fields used for training:
pred_file = os.path.join(pthmodel, model_name + "_predictors_trainidx.npz")
data_pred = np.load(pred_file)
PRED_NAMES = data_pred["PRED_NAMES"]

# Prediction grid points:
# Get MOM6 grid:
pthgrid    = config_ml["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

with xr.open_dataset(dftopo_mom) as dstopo:
  depth = dstopo['depth'].data.squeeze()

# Convert all positive values -> land (100) and ocean (<0):
HH = np.where(depth < 1.e-20, 100., -depth)

jdm, idm = HH.shape
LMsk = HH < 0

# Set mask for region and ML lat bounds:
assert HH.shape == hlat.shape, f"Shape mismatch: HH={HH.shape}, hlat={hlat.shape}"

if regn == 'north':
  DOMAIN = (hlat > lat0) & (LMsk)
elif regn == 'south':
  DOMAIN = (hlat < lat0) & (LMsk)

JG, IG = np.where(DOMAIN)

def write_nc(dfliceout, time_dnmb, A3d, regn, nctitle):
  yr1, mm1, dd1 = mtime.datevec(time_dnmb[0])[:3]
  nrecs, jdim, idim = A3d.shape

  # Days wrt to the Jan 1st:
  time_days = time_dnmb - mtime.datenum([yr1, 1, 1])
  darr_cice = xr.DataArray(A3d, dims=("time","jdim","idim"),\
                     coords={"time": time_days,\
                             "jdim": np.arange(jdim),\
                             "idim": np.arange(idim)})

  dset = xr.Dataset({"ice_thkn": darr_cice})
  dset['ice_thkn'].attrs['long_name'] = 'ice thickness'
  dset["time"].attrs = {
       "long_name": f"days since {yr1}/01/01",
       "units" : "meters",
  }

  # Add global attributes:
  dset.attrs['title']    = nctitle
  dset.attrs['source']   = "predict_gfs17_MLithkn_Ndays.py"
  dset.attrs['region']   = regn

  print(f'Dumping GLORYS interpolated ice thickness  --> {dfliceout}')
  dset.to_netcdf(dfliceout, format='NETCDF4', engine='netcdf4')


DNMB_fcst = np.arange(dnmbS, dnmbE+1)
nfcst = len(DNMB_fcst)

print(f"Start forecasts, N forecasts: {nfcst}")
for dnmb0 in DNMB_fcst:
  YR0, MM0, DD0 = mtime.datevec(dnmb0)[:3]
  rdate = YR0*10000 + MM0*100 + DD0
  print(f"Day {YR0}/{MM0}/{DD0}")

  # Derive Predictors (Input fields) 
  # hlon, hlat - grid where prediction is done (mesh025)
  # IG, JG - grid point indices where prediction is done
  # DIRS - dictionary with input / output directories
  # PREAD_NAMES - list of predictors in the right order
  # sqrt_frzdays - True: use sqrt of integrated freezing degree days
  # sst_max - max sst threshold for possible sea ice
  # standz - True: standardize predictors (False for RF)
  # regn - region of prediction
  # intgr_time - for freeze degree days, integration period, days
  # dxy - ice length scale, used for calc. ice predictors and gird point subset
  PRED, Iconc, SST = mml.construct_predictors_day(
           hlon, hlat, IG, JG, dnmb0, fyaml, PRED_NAMES,
           sqrt_frzdays, sst_max, standz, regn,
           intgr_time, Tfrz, dxy=dxy
           )
  # Construct design / predictor matrix:
  AA = np.column_stack(PRED)

  # Prediction:
  Yfcst = mlem.predict(AA)

  #Ierr = np.where(Yfcst < 0)[0]
  Yfcst[Yfcst < 0] = 0
  Yfcst[Iconc < iconc_min] = 0
  Yfcst[SST > sst_max] = 0

  if save_npz:
    pthdump = config_ml["PRED"]["pthdump"]
    os.makedirs(pthdump, exist_ok=True)
    flfcst = f"{ml_model}_ithkn_fcast_GFSv17anls_{rdate}_{regn}.npz"
    dflfcst = os.path.join(pthdump, flfcst)
    print(f"Saving fcst --> {dflfcst}")
    np.savez(dflfcst, Yfcst=Yfcst, JG=JG, IG=IG, idim=idm, jdim=jdm)

  elif save_nc:
    pthdump = config_ml["PRED"]["pthdump"]
    os.makedirs(pthdump, exist_ok=True)
    flfcst = f"{ml_model}_ithkn_GFSv17anls_mesh025_{rdate}_{regn}.nc"
    dflfcst = os.path.join(pthdump, flfcst)
    print(f"Saving fcst --> {dflfcst}")
   
    A2d = np.zeros_like(HH) * np.nan
    A2d[LMsk] = 0
    A2d[JG,IG] = Yfcst
  
    A3d = np.expand_dims(A2d, axis=0)
    time_dnmb = np.asarray([dnmb0])
    nctitle = f"ML {ml_model} predicted ithkn, GFSv17 analysis input fields"
   
    write_nc(dflfcst, time_dnmb, A3d, regn, nctitle)


if f_plt:
  import mod_colormaps as mclrmps
  from mpl_toolkits.basemap import Basemap, cm

  clrmp = mclrmps.colormap_ice_thkn()
  rmin = 0.
  rmax = 3.

  clrmp.set_bad(color=[0.1, 0.1, 0.1])
  cntr_clr = [0.9,0.,1]

  if regn == 'south':
    m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
    parallels = np.arange(-80,-10,10.)
    meridians = np.arange(-360,359.,45.)
  elif regn == 'north': 
    m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
    parallels = np.arange(50, 90, 5)
    meridians = np.arange(-360, 359., 45.)

  xh, yh = m(hlon,hlat) # GFS coords

  plt.ion()
  
  A2d = np.zeros_like(HH) * np.nan
  A2d[LMsk] = 0
  A2d[JG,IG] = Yfcst

  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.12, 0.8, 0.8])

  if regn == 'north':
    m.drawparallels(np.arange(60, 90, 10), labels=[0,0,0,0])
  elif regn == 'south':
    m.drawparallels(np.arange(-80, -50, 10), labels=[0,0,0,0])
  
  m.drawmeridians(np.arange(-180, 180, 45), labels=[0,0,0,0])
  m.drawcoastlines()
  
  img = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)

  sttl = f"ML {ml_model}, Predicted ithkn, {YR0}/{MM0:02d}/{DD0:02d}"
  ax1.set_title(sttl, fontsize=12)

  #Plot colrbar
  clb = colorbar_horiz(fig1, ax1, img, rmin=rmin, rmax=rmax, decim=2, extd='max')

  fig1.canvas.draw()

  pos_clb = clb.ax.get_position()
  bot_clb = pos_clb.y0
  pbtm = bot_clb - 0.05

  btx = f'predict_gfs17_MLithkn_Ndays.py'
  bottom_text(btx, pos=[0.02,0.02], fsz=8)


