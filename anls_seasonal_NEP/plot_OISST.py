"""
  OI SST high resolution fields
  https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html
# OpenDap to PSL data does not work on PPAN
# it works on Gaea
# I copied files PSL Linux
# [ddukhovskoy@linux256 noaa.oisst.v2.highres]$ pwd
#/Datasets/noaa.oisst.v2.highres
# ---> Niagara untrusted ---> Gaea --- gcp ---> PPAN
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import netCDF4
from netCDF4 import Dataset as ncFile
import importlib
import xarray
import yaml
from yaml import safe_load

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
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_anls_seas as manseas

YR = 1993

expt    = "seasonal_fcst"
runname = f'NEPphys_frcst_climOB_{YRS}-{MOS:02d}-e{nens:02d}'
#expt    = 'NEP_BGCphys_GOFS'
#runname = 'NEP_physics_GOFS-IC'
dnmbS   = mtime.datenum([YRS,MOS,DDS])
#dnmbR   = dnmbS + dayrun - 1

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

def read_field(furl,varnm):
  print("Reading {1} from {0}".format(furl,varnm))
  nc=ncFile(furl)
# lookup a variable
  dmm0 = nc.variables[varnm][:].data.squeeze()
  dmm = np.copy(dmm0)
  return dmm


# OpenDap to PSL data does not work on PPAN
#urlBase = 'http://psl.noaa.gov/thredds/dodsC/Datasets/noaa.oisst.v2.highres/'
#urlT    = f'sst.day.mean.{YR}.nc'
#furl = os.path.join(urlBase,urlT)
#T2d  = read_field(furl,'SST')
#ds = xarray.open_dataset(furl, chunks={})
#ds = xarray.open_dataset(furl)

pthsst = '/work/Dmitry.Dukhovskoy/data/OISST/'
flsst  = os.path.join(pthsst,f'sst.day.mean.{YR}.nc')
ds_sst = xarray.open_dataset(flsst)


