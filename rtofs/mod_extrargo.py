"""
  Extract ARGO profiles from QC output 
  log tar file
  see script 
  /scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/scripts/rtofs/get_qclog.sh
"""
import os
import numpy as np
import sys
import importlib
import datetime

def extract_argo(pthbin,flnm,Anmb=-999,Alat=-999.,Alon=-999.):
  """
    Find argo profile in the log text file
    Extract T/S
    To find 1 profile = Need to provide either argo number or lon/lat
    Otherwise - all profles are read in
  """ 

  import mod_read_ncoda as rncoda
  importlib.reload(rncoda)
  import mod_misc1 as mmisc

# ----
# Formatted text file
  fldr = pthbin+flnm

  class tsprof():
    kind = 'Argo profile'

    def __init__(self):
      self.numb  = []
      self.lon   = []
      self.lat   = []
      self.zbtm  = []
      self.ptime = []    # obs time
      self.rtime = []    # received time
      self.Tstd  = []    # T QC Std used for rejecting profile
      self.Sstd  = []    # S QC Std used for rejecting profile
      self.ZT    = []
      self.ZS    = []
      self.Tprof = []
      self.Sprof = []

    def add_data(self, numb, lon, lat, zbtm, dnmb, rnmb, Tstd, Sstd,\
                 ZT, Tdata, ZS, Sdata):
      self.numb.append(numb)
      self.lon.append(lon)
      self.lat.append(lat)
      self.zbtm.append(zbtm)
      self.ptime.append(dnmb)
      self.rtime.append(rnmb)
      self.Tstd.append(Tstd)
      self.Sstd.append(Sstd)
      self.ZT.append(ZT)
      self.ZS.append(ZS)
      self.Tprof.append(Tdata)
      self.Sprof.append(Sdata)
  
  nskip    = 0
  fid      = open(fldr,'r')
  ftmp     = False
  fsln     = False
  fhdr     = False
  ccl      = 0
  nprf     = 0
  lon      = 1.e6
  lat      = 1.e6
  ifnd     = False
  checkXY  = False
  checkN   = False
  read_all = False

  if Anmb > 0:
    checkN   = True
  elif Anmb < 0 and (abs(Alat) <= 90. and abs(Alon) <= 360.):
    checkXY  = True
  else:
    read_all = True
#    raise Exception('Either ARGO nmb or (lon,lat) should be provided')

  if checkXY:
    if Alon > 180.:
      Alon = Alon-360.

#  try:
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
      elif ifnd:
        TS.add_data(argoN, lon, lat, zbtm, dnmb, rnmb, Tstd, Sstd, \
                    ZT, Tprf, ZS, Sprf)
        if not read_all:
          break   # just 1 profile

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
        tstmp  = dmm[9]         # Argo surfacing time
        try:
          YR     = int(tstmp[0:4])
          MM     = int(tstmp[4:6])
          DD     = int(tstmp[6:8])
          HH     = int(tstmp[8:10])
          MNT    = int(tstmp[10:12])
        except:
          YR     = 9999
          MM     = 1
          DD     = 1
          HH     = 0
          MNT    = 0
        tstmpR = dmm[11]   # Argo receipt time
        try:
          YRr    = int(tstmpR[0:4])
          MMr    = int(tstmpR[4:6])
          DDr    = int(tstmpR[6:8])
          HHr    = int(tstmpR[8:10])
          MNTr   = int(tstmpR[10:12])
        except:
          YRr    = 9999
          MMr    = 1
          DDr    = 1
          HHr    = 0
          MNTr   = 0
      
        try:
          argoN  = int(dmm[-1][1:-1])
        except:
          argoN = 9999999   # number is missing 
        ptime  = datetime.datetime(YR,MM,DD,HH,MNT)
#        dnmb   = rncoda.datenum([YR,MM,DD,HH,MNT])
        dnmb   = mmisc.datenum([YR,MM,DD,HH,MNT])

        rtime  = datetime.datetime(YRr,MMr,DDr,HHr,MNTr)
#        rnmb   = rncoda.datenum([YRr,MMr,DDr,HHr,MNTr])
        rnmb   = mmisc.datenum([YRr,MMr,DDr,HHr,MNTr])

        ifnd = False
        if checkN and argoN == Anmb:
          ifnd = True
        elif checkXY and \
             np.sqrt((lat-Alat)**2+(lon-Alon)**2) < 0.01:
           ifnd = True
        elif read_all:
           ifnd = True

        if ifnd:
          print('Found profile: #{0}, lat= {1} lon={2}'.format(argoN,lat,lon))

    elif smm[:7] == "Tempera":
      try:
        Tstd = float(dmm[-1])
      except:
        Tstd = 100.  # large values *** in the file

    elif smm[:7] == "Salinit":
      try: 
        Sstd = float(dmm[-1])
      except:
        Sstd = 100.   

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
      
#  except:
#    print('Error reading file ' + fldr)

  fid.close()

  if not ifnd:
    print(' Profile {0} not found !!!'.format(Anmb))

  return TS 

