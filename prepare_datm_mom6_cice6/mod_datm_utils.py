"""
  Utility subroutines for datm file
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
#import pdb
import importlib
#import struct
import pickle
from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap
from mod_utils_fig import bottom_text


def modify_fld_nc(infile,var_name,AA):
  """
    Modify field in netCDF
    netCDF exists (template) 
  """
  ncdata = ncFile(infile,'r+')
  ncdata[var_name][:] = AA
  ncdata.close()

  return

def modify_glattr_nc(infile, var_name, var_value):
  """
    Modify global attribute in restart netCDF
    such as restart time 
  """
  ncdata = ncFile(infile,'r+')
  ncdata.setncattr(var_name, var_value)  
  ncdata.close()

  return

def addnew_glattr_nc(infile, var_name, var_value):
  """
    Add new global attribute in restart netCDF
    such as information about restart files/ authors, etc/

  Note: the code freezes up trying to add a new attribute
  Need to debug why
  """
  ncdata = ncFile(infile,'r+')
  ncdata.setncattr(var_name, var_value)  
  ncdata.close()

  return



def read_ncfile(flname, fldname, fsilent=False):
  if not fsilent:
    print('reading ' + 'fldname' + ' from ' + flname)
  ncdata = ncFile(flname, 'r')
  AA    = ncdata[fldname][:].data.squeeze()
  return AA 

def add_fld2D_ncfile(flname, fldname):
  """
    Add new variable to existing netCDF
  """
  print('Adding ' + fldname + ' --> ' + flname)
  infile = ncFile(flname, 'r+')
  nx_nc  = infile.dimensions['ni'].name
  ny_nc  = infile.dimensions['nj'].name
  newvar = infile.createVariable('coszen','f8',(ny_nc,nx_nc))
  infile.close() 

  return
 
def convert_sec2date(nsecs, ref_datenum):
  """
  Convert seconds since ref_datenum to datenum
  Return: list of [dnmb, HR, YY, MM, DD, HR, minutes]
  """
  import mod_time as mtime
  Ndays     = nsecs/86400.
  time_dnmb = ref_datenum + Ndays
  DV        = mtime.datevec(time_dnmb)
  YR        = DV[0]
  MM        = DV[1]
  DD        = DV[2]
  HR        = DV[3]
  MN        = DV[4]

  date_vec = [time_dnmb, YR, MM, DD, HR, MN]

  return date_vec

def convert_dnmb2sec(dnmb, ref_dnmb):
  """
  Compute N seconds from reference date ref_dnmb
  """
  Ndays = dnmb-ref_dnmb
  Nsecs = Ndays*86400.

  return Nsecs

def datm_newfile(fl_tmplt, fl_new, fovrd=False):
  """
    Create a new netCDF file from some template file
    fl_tmplt - template, does not change this file
    fl_new   - new file copied from the template for editing
    if fovrd - rewrite existing fl_new
    otherwise - do not do anything
  """
  import shutil

  if os.path.exists(fl_new):
    print(fl_new + ' exists, will modify this file')
    if fovrd:
      print(fl_new + ' will be overode')
    else:
      return 

  print(fl_tmplt + ' ---> ' + fl_new)
  shutil.copy(fl_tmplt, fl_new)

  return


