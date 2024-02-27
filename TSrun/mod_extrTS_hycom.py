"""
  Extract T/S profiles from HYCOM/ROFS archv nonformatted binary archv
"""
import os
import numpy as np
import sys
import importlib
import datetime
import mod_utils as mutil
import mod_read_hycom as mhycom

def extract_1prof(pthbin, flnm, fld, lon0, lat0, LON, LAT,\
                  fthknss=True, KDM=-1):
  """
    Find profile closest to the location (lon0,lat0)
    fld - name of the field to be extracted, should match *.b 
    LON, LAT - RTOFS/HYCOM grid coordinates
    fthknss - if true, extract layer thiknesses and ZZ, ZM profiles
            corresponding to fld
  """
  ii0, jj0 = mutil.find_indx_lonlat(lon0,lat0,LON,LAT)

  fina = pthbin+flnm+'.a'
  finb = pthbin+flnm+'.b'
  rg   = 9806.
  huge = 1.e20

  if KDM <= 0:
    IDM, JDM, KDM = mhycom.hycom_dim(fina,finb)
  else:
    JDM = LON.shape[0]
    IDM = LON.shape[1]

  ZZ = []
  ZM = []
  dH = []
  if fthknss:
    for kk in range (1,KDM+1):
      F,nn,mm,ll = mhycom.read_hycom(fina,finb,'thknss',rLayer=kk,finfo=False)
      F = F/rg
      F[np.where(F>huge)] = np.nan
      F[np.where(F<0.001)] = 0.
      dH.append(F[jj0,ii0])

    dH = np.array(dH)
    ZZ, ZM  = mhycom.zz_zm_fromDP1D(dH,finfo=False)

  F1D = []
  for kk in range (1,KDM+1):
    F,nn,mm,ll = mhycom.read_hycom(fina,finb,fld,rLayer=kk,finfo=False)
    F[np.where(F>huge)] = np.nan
    F1D.append(F[jj0,ii0])
    
  F1D = np.array(F1D)

  return F1D, ZZ, ZM, dH

