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


iconc_new = 0.3016
iconc_old = 0.0
ithk_new  = 2214.57    # kg m-2
ithk_old  = 2327.25

# Find change in ice = dlt_iconc*ithk_old + dlt_ithk*iconc_new
# same for snow 
dlt_iconc = iconc_new - iconc_old
dlt_ithk  = ithk_new - ithk_old
dlt_ice   = dlt_iconc * ithk_old + dlt_ithk * iconc_new

# Better way: compute mass*fract_area 
dlt_ice   = (ithk_new*iconc_new - ithk_old*iconc_old)
# Using bluk ice S:
ice_salin = 3.35
dlt_salt = dlt_ice * ice_salin  # kg m-2 * g kg-1

# Enthalpy for BL99 thermodynamics - Turner & Hunke, 2015:
# The enthalpy is the energy required to melt the ice completely and raise its 
# temperature to 0C.
# From /gpfs/f5/cefi/scratch/Dmitry.Dukhovskoy/src/FMS/constants4
# gfdl_constants.fh
c0 = 2.1060e3 # specific heat of fresh ice 0C, J/[kg *degC]  [Q C-1 ~> J kg-1 degC-1]
L0 = 3.34e5 # latent heat of fusion of fresh ice at 0C [Q ~> J kg-1]
Tm = -0.1 # ice melting T - function of ice (S)
Cw = 3989.24495292815 # specific heat capacity of brine or sea water [Q C-1 ~> J kg-1 degC-1]
rho_ice = 990.
Tice = -1.8
dlt_ice = 0.6482e3  # kg/m2 - ice change


# Enthalpy of ice, J/m3
q_Jm3 = -rho_ice*(c0*(Tm-Tice) + L0*(1.-Tm/Tice) - Cw*Tm)

# Enthalpy of ice, J/kg
q_Jkg = -(c0*(Tm-Tice) + L0*(1.-Tm/Tice) - Cw*Tm)

# dlt heat:
dlt_heat = q_Jkg * dlt_ice
print(f'enthalpy J/m3 = {q_Jm3:12.1f}, J/kg = {q_Jkg:12.1f}')
print(f'dlt_ice kg/m2 = {dlt_ice:12.1f}, dlt_heat J/m2 = {dlt_heat:12.1f}')
 
