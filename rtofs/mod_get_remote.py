"""
  Read remote data from hycom.org
"""
import os
import numpy as np
import sys
import netCDF4 as nc
import importlib
from netCDF4 import Dataset as ncFile

def read_gofs31_ts3z(YR,MM,DD):
  """
    GOFS 3.1: 41-layer HYCOM + NCODA Global 1/12° Analysis
    OPENDAP: //tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/ts3z
    NetcdfSubset: //ncss.hycom.org/thredds/ncss/grid/GLBy0.08/expt_93.0/ts3z
    WMS: //wms.hycom.org/thredds/wms/GLBy0.08/expt_93.0/ts3z
    WCS: //wcs.hycom.org/thredds/wcs/GLBy0.08/expt_93.0/ts3z
  """
  rdate = '{0}{1:02d}{2:02d}12'.format(YR,MM,DD)
  url_data = 'http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/ts3z/'
#  url_data = 'https://ncss.hycom.org/thredds/ncss/grid/GLBy0.08/expt_93.0/ts3z/' 
  url_data = 'http://data.hycom.org/datasets/GLBy0.08/expt_93.0/data/hindcasts/{0}/'.\
             format(YR)
  file_url = url_data + 'hycom_glby_930_'+rdate+'_t000_ts3z.nc'
  DA = ncFile(file_url)


