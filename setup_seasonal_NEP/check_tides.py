"""
  Check tides imposed at the OB
  OBC_TIDE_CONSTITUENTS = "M2,S2,N2,K2,K1,O1,P1,Q1,MM,MF"
  MF - The near-fortnightly tide Mf, of period 13.66 d, is the largest of  
       the zonally symmetric, long-period tides. Like all the long-period lunar tides, 
       it may be thought of as a time-varying modulation of the Earth permanent tide M0.
       The MF tide is generally 2 cm or less.
  
"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
#import pickle
import matplotlib.pyplot as plt
from yaml import safe_load

import mod_utils_ob as mutob
importlib.reload(mutob)


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
sys.path.append('./seasonal-workflow')
from boundary import Segment
import mod_time as mtime
import mod_mom6 as mmom6
import mod_utils as mutil
from mod_utils_fig import bottom_text

varnm  = 'ssh'

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

fconfig = 'config_nep.yaml'
with open(fconfig) as ff:
  config = safe_load(ff)

# MOM6 NEP topo/grid:
run_name   = 'seasonal_fcst_daily'
pthtopo    = gridfls['MOM6_NEP'][run_name]['pthgrid']
fgrid      = gridfls['MOM6_NEP'][run_name]['fgrid']
ftopo_mom  = gridfls["MOM6_NEP"][run_name]["ftopo"]
outdir     = gridfls['MOM6_NEP'][run_name]['pthoutp']
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
dftopo_mom  = os.path.join(pthtopo, ftopo_mom)
LONM, LATM  = mmom6.read_mom6grid(dfgrid_mom)
HHM         = mmom6.read_mom6depth(dftopo_mom)

segments = [ Segment(1, 'north', hgrid, output_dir=outdir),
             Segment(2, 'east',  hgrid, output_dir=outdir),
             Segment(3, 'south', hgrid, output_dir=outdir),
             Segment(4, 'west',  hgrid, output_dir=outdir)]

nOB = len(segments)

pthtides='/work/Dmitry.Dukhovskoy/NEP_input/tides_OBs/'


nsgm = 1
const=9
flnm = f'tz_{nsgm:03d}.nc'
#flnm = f'tz_{nsgm:03d}_test.nc'
dfinput = os.path.join(pthtides, flnm)
dset = xarray.open_dataset(dfinput)
AA = dset[f'zamp_segment_{nsgm:03d}'].isel(time=0, constituent=const-1).data.squeeze()



plt.ion()
fig1 = plt.figure(1,figsize=(9,8))

fig1.clf()
ax1  = plt.axes([0.1, 0.3, 0.8, 0.6])
ax1.plot(AA)

sttl = f'Tidal amplitude, const={const} at the OB segment {nsgm}'
ax1.grid('on')
ax1.set_title(sttl)
ax1.set_xlabel('Dist., km')

btx='check_ssh_dailyOB.py'
bottom_text(btx)



