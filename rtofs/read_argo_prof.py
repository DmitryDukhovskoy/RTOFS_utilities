"""
  Read ARGO profiles from QC output 
  log tar file
  see script 
  /scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/scripts/rtofs/get_qclog.sh
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
import datetime
#import pickle
#from netCDF4 import Dataset as ncFile
import matplotlib.colors as colors
import matplotlib.mlab as mlab
#from matplotlib.patches import Polygon
#from matplotlib.colors import ListedColormap
#from mpl_toolkits.basemap import Basemap, shiftgrid

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')

import mod_read_ncoda as rncoda
importlib.reload(rncoda)

rdate   = '20230123'
pthbin  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/rtofs.'+\
          rdate+'/ocnqc_logs/profile_qc/'

yr, mo, mday, hr = rncoda.parse_rdate(rdate)

time_n00  = datetime.datetime(yr,mo,mday,hr,0,0)
time_n24  = time_n00 - datetime.timedelta(days=1)
rdate_n24 = time_n24.strftime('%Y%m%d%H') 


# ----
# Formatted text file
flnm = 'prof_argo_rpt.'+rdate_n24+'.txt' 
fldr = pthbin+flnm


class tsprof():
  kind = 'Argo profile'

  def __init__(self):
    self.numb  = []
    self.lon   = []
    self.lat   = []
    self.zbtm  = []
    self.ptime = []
    self.ZT    = []
    self.ZS    = []
    self.Tprof = []
    self.Sprof = []

  def add_data(self, numb, lon, lat, zbtm, dnmb, ZT, Tdata, ZS, Sdata):
    self.numb.append(numb)
    self.lon.append(lon)
    self.lat.append(lat)
    self.zbtm.append(zbtm)
    self.ptime.append(dnmb)
    self.ZT.append(ZT)
    self.ZS.append(ZS)
    self.Tprof.append(Tdata)
    self.Sprof.append(Sdata)
  
nskip = 0
fid = open(fldr,'r')
ftmp = False
fsln = False
fhdr = False
ccl  = 0
nprf = 0
lon  = 1.e6
lat  = 1.e6
try:
# Find length of the input file:
  fid.seek(0,2)
  fend = fid.tell()
# Read header then data
  fid.seek(0)
  fpos = 1

  while fpos < fend:
    dmm = fid.readline().split()
    if len(dmm)==0: 
      continue
    elif len(dmm[0].strip()) == 0:
      print('Empty string, skipping ...')
      continue
  
    smm = str(dmm[0])
    fpos = fid.tell()
    if fpos == fend:
      break

    if smm[:5] == "-----":
      fhdr = True

    elif smm[:7] == "Profile":
# Save profiles
      if nprf == 0:
        TS = tsprof()
      else:
        TS.add_data(argoN, lon, lat, zbtm, dnmb, ZT, Tprf, ZS, Sprf)

      nprf += 1
      print('Reading profile # {0}'.format(nprf))
      ZT    = []
      ZS    = []
      Tprf  = []
      Sprf  = []
      if dmm[2][:3] == "Lat":
        lat    = float(dmm[3])
        lon    = float(dmm[5])
        zbtm   = float(dmm[7])
        tstmp  = dmm[9]
        YR     = int(tstmp[0:4])
        MM     = int(tstmp[4:6])
        DD     = int(tstmp[6:8])
        HH     = int(tstmp[8:10])
        MNT    = int(tstmp[10:12])
        try:
          argoN  = int(dmm[-1][1:-1])
        except:
          argoN = 9999999   # number is missing 
        ptime  = datetime.datetime(YR,MM,DD,HH,MNT)
        dnmb   = rncoda.datenum([YR,MM,DD,HH,MNT])
    elif dmm[0] == "Depth":
      fhdr = False
      if dmm[1] == "Temp":
        ftmp = True
        fsln = False
      else:
        ftmp = False
        fsln = True
      continue

# Reading data:      
    if not fhdr:
      zz  = -1.*float(dmm[0])
      val = float(dmm[1]) 

      if ftmp:
        ZT.append(zz)
        Tprf.append(val)
      else:
        ZS.append(zz)
        Sprf.append(val)
    
except:
  print('Error reading file ' + fldr)

fid.close()


