"""
  Final step in preparing MOM ocean restart
  Combin all fields and write a netcdf file
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import datetime
import pickle
import yaml
from netCDF4 import Dataset as ncFile
import xarray as xr

#PPTHN = '/home/Dmitry.Dukhovskoy/python'
PPTHN = []
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_read_hycom as mhycom
import mod_colormaps as mcmp
import mod_mom6 as mom6util
import mod_regmom as mrgm
#import mod_valid_utils as mvutil
importlib.reload(mcmp)

nrun       = "GOFS3.0"
nexpt      = 19.0
expt       = f"{nexpt:2.1f}"
YR         = 1993
jday       = 1
hr         = 15
grid_shape = 'symmetr'   # MOM grid: symmetr/nonsymmetr

# Restart date:
YR_r    = 1993
jday_r  = 1
hr_r    = 0


pthout  = '/work/Dmitry.Dukhovskoy/data/mom6_nep_restart/'
pthtmp  = '/work/Dmitry.Dukhovskoy/data/mom6_nep_tmp/'
pthgofs = '/work/Dmitry.Dukhovskoy/data/GOFS3.0/expt_19.0/'

dnmb   = mtime.jday2dnmb(YR, jday)
dv     = mtime.datevec(dnmb)
ds     = mtime.datestr(dnmb)
MM     = dv[1]
DM     = dv[2]
rdate  = f"{YR}{MM:02d}{DM:02d}"

dnmb_r = mtime.jday2dnmb(YR_r, jday_r)
dv_r   = mtime.datevec(dnmb_r)
ds_r   = mtime.datestr(dnmb_r)
MM_r   = dv_r[1]
DM_r   = dv_r[2]
rdate_r= f"{YR_r}{MM_r:02d}{DM_r:02d}{hr_r:02d}"

flnm_rst  = f"MOM.res.{rdate_r}.nc"
dflnm_rst = os.path.join(pthout, flnm_rst) 

print(f" output MOM restart file: {dflnm_rst}")

# Read MOM6 NEP domain nx=342, ny=816
with open('pypaths_gfdlpub.yaml') as ff:
  dct = yaml.safe_load(ff)

pthgrid_mom = dct["MOM6_NEP"]["test"]["pthgrid"]
ftopo_mom   = dct["MOM6_NEP"]["test"]["ftopo"]
fgrid_mom   = dct["MOM6_NEP"]["test"]["fgrid"]
dftopo_mom  = os.path.join(pthgrid_mom, ftopo_mom)
dfgrid_mom  = os.path.join(pthgrid_mom, fgrid_mom)


def load_interp(pthtmp, fldiout, grid_shape, rdate):
  flintrp    = f"gofs2mom_nep_vrtz-{fldiout}_{grid_shape}_{rdate}.pkl"
  dflintrp   = os.path.join(pthtmp,flintrp)
  print(f"Loading {dflintrp}")
  with open(dflintrp, 'rb') as fid:
    Ai = pickle.load(fid)

  return Ai

lonh, lath  = mom6util.read_mom6grid(dfgrid_mom, grid=grid_shape, grdpnt='hgrid')
lonq, _     = mom6util.read_mom6grid(dfgrid_mom, grid=grid_shape, grdpnt='ugrid') 
_, latq     = mom6util.read_mom6grid(dfgrid_mom, grid=grid_shape, grdpnt='vgrid')

time_var = np.atleast_1d(dnmb)
nvlrs    = 75
lonh_var = lonh[0,:].squeeze()
lath_var = lath[:,0].squeeze()
lonq_var = lonq[0,:].squeeze()
latq_var = latq[:,0].squeeze()
lrs_var  = [x*1. for x in range(0,nvlrs)]


Ai = load_interp(pthtmp, 'uvel', grid_shape, rdate)
Ai = np.expand_dims(Ai, axis=0)
coords = {
  "Time": time_var,
  "Layer": lrs_var,
  "lath": lath_var,
  "lonq": lonq_var
}
udim = ['Time','Layer','lath','lonq'] 
uattrs = {'long_name': 'Zonal velocity',
         'units': 'm s-1'}
udarray = xr.DataArray(Ai, coords, dims=udim, attrs=uattrs, name="v")

Ai = load_interp(pthtmp, 'vvel', grid_shape, rdate)
Ai = np.expand_dims(Ai, axis=0)
coords = {
  "Time": time_var,
  "Layer": lrs_var,
  "latq": latq_var,
  "lonh": lonh_var
}
vdim = ['Time','Layer','latq','lonh'] 
vattrs = {'long_name': 'Meridional velocity',
         'units': 'm s-1'}
vdarray = xr.DataArray(Ai, coords, dims=vdim, attrs=vattrs, name="u")

Ai = load_interp(pthtmp, 'salin', grid_shape, rdate)
Ai = np.expand_dims(Ai, axis=0)
coords = {
  "Time": time_var,
  "Layer": lrs_var,
  "lath": lath_var,
  "lonh": lonh_var
}
sdim = ['Time','Layer','lath','lonh']
sattrs = {'long_name': 'Salinity',
         'units': 'PPT'}
sdarray = xr.DataArray(Ai, coords, dims=sdim, attrs=sattrs, name="Salt")

Ai = load_interp(pthtmp, 'temp', grid_shape, rdate)
Ai = np.expand_dims(Ai, axis=0)
coords = {
  "Time": time_var,
  "Layer": lrs_var,
  "lath": lath_var,
  "lonh": lonh_var
}
tdim = ['Time','Layer','lath','lonh']
tattrs = {'long_name': 'Potential Temperature',
         'units': 'degC'}
tdarray = xr.DataArray(Ai, coords, dims=tdim, attrs=tattrs, name="Temp")

Ai = load_interp(pthtmp, 'lrthknss', grid_shape, rdate)
Ai = np.expand_dims(Ai, axis=0)
coords = {
  "Time": time_var,
  "Layer": lrs_var,
  "lath": lath_var,
  "lonh": lonh_var
}
ldim = ['Time','Layer','lath','lonh']
lattrs = {'long_name': 'Layer Thickness',
         'units': 'm'}
ldarray = xr.DataArray(Ai, coords, dims=ldim, attrs=lattrs, name="h")


DS = xr.merge([udarray, vdarray, ldarray, tdarray, sdarray])
DS.Time.attrs={'long_name': 'Time',
               'cartesian_axis': "T",
               'units': "days"}
DS.Layer.attrs={'long_name': "Layer z-rho",
                'cartesian_axis': "Z",
                'units': "meter",
                'positive': "up"}
DS.lath.attrs={'long_name':  "Latitude",
               'cartesian_axis': "Y",
               'units': "degrees_north"}
DS.latq.attrs={'long_name':  "Latitude",
               'cartesian_axis': "Y",
               'units': "degrees_north"}
DS.lonh.attrs={'long_name':  "Longitude",
               'cartesian_axis': "X",
               'units': "degrees_east"}
DS.lonq.attrs={'long_name':  "Longitude",
               'cartesian_axis': "X",
               'units': "degrees_east"}
DS.attrs={'history':f'Created from GOFS3.0 {rdate}'}
DS.info()

print(f"\nWriting netcdf ---> {dflnm_rst}")
DS.to_netcdf(dflnm_rst, mode='w', unlimited_dims='Time')


f_chck = False
if f_chck:
  plt.ion()

  ilr = 1
  aa = Ai[ilr-1,:,:].squeeze()
  stl = f"Lr = {ilr}"
  idm = aa.shape[1]
  jdm = aa.shape[0]
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  im1 = ax1.pcolormesh(aa, \
                   vmin=10., \
                   vmax=33.)

  ax1.axis('scaled')
  ax1.set_xlim([0, idm-1])
  ax1.set_ylim([0, jdm-1])
  ax1.set_title(stl)

  ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                     0.02, ax1.get_position().height])
  clb = plt.colorbar(im1, cax=ax2, orientation='vertical', extend='both')





