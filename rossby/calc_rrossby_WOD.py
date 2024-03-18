"""
  1st baroclinic Rossby radius 
  Numerically solve eig/value problem
  Use WOD T/S profiles
 
  Use bottom topography for depths deeper than WOA last depth level (-5500m)
   Topography interpolated from ETOPO2
 
 
  Calculate N2 following Chelton 1996
  N2 is calculated at the middle of the depth interval
  density is adiabatically adjusted to the mid-grid depth

  extract profiles and calc R. radius & Phase speed of the gravity wave
  for the 1st baroclinic mode
  for time period YR1-YR2

  The radius and phase speed are saved in a binary files
  at /data/ncei2/OCL/sys.inf/
  from all probe types
  Position in the files: The file is indexed by unique cast number 
  so position NNN in the file will contain the calculated value 
  for unique cast NNN 

  Dmitry Dukhovskoy, NOAA OAR PSL
  March 2024

"""
import os
import numpy as np
import matplotlib.pyplot as plt
#import torch
import sys
import pdb
import importlib
#import netCDF4
#from netCDF4 import Dataset as ncFile
import timeit
from os.path import exists
import yaml
from netCDF4 import Dataset as ncFile

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
import mod_misc1 as mmsc1


# -------------------------------------------
#
#   USER OPTIONS: select time, flags 
#   
#   Flags: 
#   f_save  = 1 - saves Rossby R and C phase to sys.inf/. binary
#                  for computing climatologies 
#                 and saves list of processed observations
#                 The list of processed observations is saved by years
#                 in separate files in a scratch directory 
#                 These are temporary files that can be deleted later
#           = 0 - no output saved, use for test/debug
#
#   f_load  = 1  Load saved processed profiles info
#                i.e. lists of processed sequential records
#                (for all probes) are loaded and 
#                records to be processed are updated
#                so that processed records are not processed again
#                if file does not exist, f_load is turned off
#           = 0  All observations in all probes will be processed
#                from the beginning of the list of selected
#                observations that passed the QC flag for
#                computing the Radius
#
#  f_addWOA = 1  Use WOA18 0.25 grid T/S profiles
#                to fill missing portions of observed T/S profiles 
#                in the deep layers
#           = 0  Use only observed T/S profiles with possible
#                gaps in the deep layers 
#                Incomplete profiles may result in big errors in 
#                Rossby R estimates
#                When WOA is not used - more strict QC is 
#                used which drastically reduces
#                the number of WOD profiles 
# -------------------------------------------
YR1 = 2002
YR2 = 2002
f_save   = 1    # = 1 - saves Rossby R and C phase to sys.inf/. binary
f_load   = 1    # Load saved processed profiles info
f_addWOA = 0

pthout = '/data/ncei2/work/dmitry.dukhovskoy/data_output/'
pthrsb = '/data/ncei2/OCL/sys.inf/'  # dir for Rossby r. output
btx = 'calc_rrossby_WOD.py'

if f_addWOA > 0:
  sfx = '_woa'
else:
  sfx = ""

# Keep info of processed profiles
# for each Probe
"""
class ProfProcessed:
  def __init__(self,recnmb=-999,dnmb=-999,lon=-999.,lat=-999.,\
               crsnmb=-999.,uniqnmb=-999.):
    self.recnmb[0] = recnmb  # seq record number
    self.dnmb[0]   = dnmb
    self.lon[0]    = lon
    self.lat[0]    = lat
    self.crsnmb[0] = crsnmb
    self.uniqnmb[0]= uniqnmb   # unique profile number
"""

class ListRcrds:
# List of sequential records of the casts for each probe
  def __init__(self,ListRcrd):
    self.list_rcrd = ListRcrd

import mod_utils_wod as uwod
importlib.reload(uwod)
import mod_utils_rossby as ursb
importlib.reload(ursb)
import mod_calc_rrosby as clcr
importlib.reload(clcr)
import mod_solver as uslv
importlib.reload(uslv)
import mod_plot_rsb as uplt
importlib.reload(uplt) 

# Read WOA18 grid
# lon_woa/lat_woa are identical to
# LONW, LATW from interpolated TOPO
if f_addWOA == 1:
  lon_woa, lat_woa, Z_woa = ursb.read_WOAgrid025()

RCRDS = {}   # dict with the list of seq. records for all probes
PRBS=["CTD","OSD","MRB","PFL","DRB"]
#PRBS = ["DRB","MRB"]
print(' Time range: {0}-{1} \n Probes to process :'.format(YR1,YR2))
for prbnm in PRBS:
  print('====> '+prbnm)

if f_addWOA == 1:
  print(' WOA18 climatology used to fill gaps in deep layers')
else:
  print(' WOA18 climatology NOT used to fill gaps')

for YR in range(YR1,YR2+1):
# Processed profiles
  fout = pthout+'processed_casts{0}{1}.pkl'.format(YR,sfx)    
  fprbs = pthout+'list_seqrcrds{0}.pkl'.format(YR)

# List of all sequential records per probe 
# That fall within the time range
# The list takes time to generate
# for some probes (e.g., CTD), thus
# save it once created and load if needed
  if exists(fprbs):
    print('Load List of seq.records '+fprbs)
    with open(fprbs,'rb') as fid:
      RCRDS = pickle.load(fid)
  else:
    for prbnm in PRBS:
      print('Finding casts for {0} {1}'.format(prbnm,YR))
      Lst_rcrd = ursb.find_recs_prcs(prbnm,YR1=YR1,YR2=YR2)
      RCRDS[prbnm] = Lst_rcrd   # Save seq.records to process by probe 
    
    print('Search of seq.records done, dumping --> '+fprbs)
    with open(fprbs,'wb') as fid:
      pickle.dump(RCRDS,fid)

# For statistics and process control, keep initial # of records
  for prbnm in PRBS:
    Lst_rcrd = RCRDS.get(prbnm)
    nprbnm = 'n'+prbnm
    RCRDS[nprbnm] = len(Lst_rcrd) # initial # of records to process

# Processed records:
  if f_load > 0:
# Read in saved processed file
    if exists(fout):
      print('Loading saved processed records #'+fout)
      with open(fout,'rb') as fid:
        PRCSD = pickle.load(fid)
    else:
      print('!!!!! Processed records not saved '+fout)
      print('!!!!!    Initializing PRCSD  !!!!!!!\n')
      f_load = 0

  if f_load == 0:
# Initiate list of processed sequential casts by probe
    PRCSD={}
    for prbnm in PRBS:
      PRCSD[prbnm]=[]
      qcnm = prbnm+'QC'
      PRCSD[qcnm]=[]
#
# Next: update list of records to process excluding processed:
  if f_load > 0:
    RCRDS = ursb.update_rcrds(PRCSD,RCRDS)

  for prbnm in PRBS:
    nprbnm = 'n'+prbnm
    Nrc0   = RCRDS.get(nprbnm)
    Nprcsd = len(PRCSD.get(prbnm))
    print('===> {0} fraction to process: {1:.2f}%'.\
          format(prbnm,(1.-float(Nprcsd)/float(Nrc0))*100.))
  

# ------------
# Read interpolated ETOPO2:
# -------------
HW,LONW,LATW = ursb.read_topo_WOA025()

# ==========================================
#
#        START
#
# ==========================================
print('\n  Starting the calculation of the Rossby Radius \n')
for prbnm in PRBS:
  print('  Probe = '+prbnm)
  Lst_rcrd = RCRDS.get(prbnm) 
# Initial record length:
  nprbnm = 'n'+prbnm
  Nrc0   = RCRDS.get(nprbnm)  # initial # of records
  Nrc    = len(Lst_rcrd)      # # of records left
  Nprcsd = len(PRCSD.get(prbnm))

  for irc in range(Nrc):
# Read T/S profiles and obs. depths
    krcrd = Lst_rcrd[irc]
    done_tot = float(Nprcsd)/float(Nrc0)+float(irc)/float(Nrc0)
    print('Reading {0}  seq.rcrd# {1}, done: {2:.3f}%, overall {3:.3f}%'.\
          format(prbnm,krcrd,float(irc)/float(Nrc)*100.,done_tot*100.))

    PRBINFO, PRB = uwod.read_header(prbnm,krcrd)
    Zobs = uwod.read_obsdepths(PRBINFO,PRB,fecho=0)   # observed depth levels
    Tobs = uwod.read_varobs(PRBINFO,PRB,'Temperature',fecho=0) # at obs. depths
    Sobs = uwod.read_varobs(PRBINFO,PRB,'Salinity',fecho=0) # at obs. depths
    Tobs[np.where(np.abs(Tobs)>1.e4)]=np.nan
    Sobs[np.where(np.abs(Sobs)>1.e4)]=np.nan
# Error T/S values:
    Tobs[np.where(Tobs>40.)]   = np.nan
    Tobs[np.where(Tobs<-2.0)]  = np.nan
    Sobs[np.where(Sobs>50.)]   = np.nan
    Sobs[np.where(Sobs<1.e-6)] = np.nan
#
    if np.min(Zobs) >= 0.:
      Zobs=-Zobs
#
# Find local depth
    lon0 = PRB.lon
    lat0 = PRB.lat
    zbtm = ursb.find_depth(HW,LONW,LATW,lon0,lat0)

# When observed profile is too short wrt the local depth
# this results in errors in the Radius calculation - it is bigger
# than when computed with the full profiles 
# (when deep weakly stratified part is missing)
# one option is for incomplete profiles to add the WOA18 for missing part
    if f_addWOA == 1:
# Run QC before adding WOA with addWOA flag
# then after without addWOA flag
      f_QC = ursb.run_QC(Zobs,Tobs,Sobs,zbtm,PRB, f_addWOA=f_addWOA)
# Update processed list if QC failed
      if f_QC > 0:
        PRCSD = ursb.update_PRCSD(prbnm,f_QC,krcrd,PRCSD)
        print('<<>>> Skipping seqrcrd {0} {1} QC flag={2} zbtm={3:.1f}m'.\
            format(prbnm,krcrd,f_QC,zbtm))
        continue

      Tin = Tobs.copy()
      Sin = Sobs.copy()
      Zin = Zobs.copy()
      Tobs, Sobs, Zobs = ursb.add_WOA(Tobs,Sobs,Zobs,PRB,\
                                      LONW,LATW,Z_woa)

# In some cases, WOA bottom depth can be shallower
# or onland
# than use the deepest obs. depth as bottom with small offset
    if zbtm > Zobs[-1]:
      zbtm = Zobs[-1]-1.5

# QC check of T/S profiles
    f_QC = ursb.run_QC(Zobs,Tobs,Sobs,zbtm,PRB)

# Update processed list if QC failed
    if f_QC > 0:
      PRCSD = ursb.update_PRCSD(prbnm,f_QC,krcrd,PRCSD)
      print('>>>> Skipping seqrcrd {0} {1} QC flag={2} zbtm={3:.1f}m'.\
            format(prbnm,krcrd,f_QC,zbtm))
      continue

# Vertical interp of T/S into Z depths
# Main purpose to fill missing observed values in both profiles
# at observed depths (default) 
# Other options: extend profiles to the surface (add_srf=1)
# interpolate onto a high-res. vertical grid for better 
# matrix inversion when solving the Sturm-Liouville problem (fintrp=1)
# Returned:
# ZZ = depths of the grid cell interfaces, i.e. depths corresponding
# to Tzz, Szz values
# If fintrp = 0 (or not specified), ZZ = Zobs 
# 
#    ZZ = ZZ
#   Tzz, Szz, ZZ = ursb.interp_gaps_TS(Zobs,Tobs,Sobs,zbtm,add_srf=0)
    Tzz, Szz, ZZ = ursb.interp_TS_zlev(Zobs,Tobs,Sobs,zbtm)

    f_plt = 0
    if f_plt:
      uniqnmb = PRB.uniqnmb
      ctl = '{1}, U#{2} Zbtm={0:.1f}'.format(zbtm,prbnm,uniqnmb)
#      uplt.plot_prof(Zobs,Tobs,ctl=ctl) # original
      ctl2='T Interpolated'
      uplt.plot_2prof(Zobs,Tobs,ZZ,Tzz,ctl1=ctl,ctl2=ctl2,fgn=1) 
      bottom_text(btx)

      ctl3 = '{1}, U#{2} Zbtm={0:.1f}'.format(zbtm,prbnm,uniqnmb)
      ctl4='S Interpolated'
      uplt.plot_2prof(Zobs,Sobs,ZZ,Szz,ctl1=ctl3,ctl2=ctl4,fgn=2) 
      bottom_text(btx)

# Compute N2 profile:
# Neutral method, centered-difference - compute in middle grid points
# Z_phi - depths where N2 is calculated
# The N2 profiles in most cases are not full water-columns profiles
# but only over the section of the water column that was observed
# Filling down to the bottom is performed when e/value problem is solved
    N2, Z_phi = clcr.calcN2_1D(Tzz,Szz,ZZ,lat0,info=0)

# Fill the gap from the deepest obs to the bottom
    N20 = N2.copy()
    Z_phi0 = Z_phi.copy()
    ZZ0 = ZZ.copy()
    N2, Z_phi, ZZ = clcr.intrpN2_bottom(N2,Z_phi,ZZ,zbtm)

# Compute Rossby R
# for 1 profile
# Rrsb - 1st barocl Rossby radius, km
# Cphase - 1st baroclinic gravity wave phase speed m/s
# Phi - eigenvector corresponding to the 1st e/value
    Rrsb, Cphase, Phi = clcr.solve_SturmLiouv_1D(N2,Z_phi,ZZ,zbtm,lat0)
    print('---> lon={0:.1f}E, lat={1:.1f}N, Zbtm={2:.1f}'.\
           format(lon0,lat0,zbtm))
    print('---> Rossby R={0:.1f} km, Cphase = {1:.1f} m/s\n'.format(Rrsb,Cphase))

# Test - plot N profiles
# Convert to cycles/hr
    f_testN = 0
    if f_testN == 1:
      uplt.plot_checkN2Phi(N2,ZZ,Z_phi,Phi,PRB,zbtm,prbnm,krcrd,Rrsb,Cphase)
      bottom_text(btx, pos=[0.01, 0.02])

# Update processed list
    PRCSD = ursb.update_PRCSD(prbnm,f_QC,krcrd,PRCSD)

# Write Radius 
    uniqnmb = PRB.uniqnmb
    if f_save == 1:
      frsb = pthrsb+'RossbyRadius'+sfx+'.s'
      val = Rrsb.item()
      print(' Writing U#{0} Rrsb--> {1}'.format(uniqnmb,frsb))
      ursb.write_1output(frsb,val,PRB)

      fcphase = pthrsb+'Cphase'+sfx+'.s'
      val = Cphase.item()
      print(' Writing Cphase--> '+fcphase)
      ursb.write_1output(fcphase,val,PRB)
 
# Save processed obs. casts only if binary output is saved:
      with open(fout,'wb') as fid:
        print(' Saving processed recrods # in '+fout)
        pickle.dump(PRCSD,fid)
     
    print('=================\n')    

if f_save == 1:
  print(' Final: Processed records ')
  with open(fout,'wb') as fid:
    print(' Saving processed recrods # in '+fout)
    pickle.dump(PRCSD,fid)

print('Done ')

