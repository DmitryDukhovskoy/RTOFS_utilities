"""
  Plot NCODA increment input fields
  in layers - new ncoda approach
  used in ncoda_archv_lyrinc.f
  old:
  grdlvl_lyr*datafld  - mid-points of the forecast layer pressures in layerspace(Pa)
  New:
  *_0024_fcstfld - background/f/cast fields from HYCOM
  *_0000_analinc - NCODA increments in layers
  *_0000_analfld - NCODA updated (background+increment) in layers

  lyrprs_lyr*_0024_fcstfld - mid-points of the bacground (f/cast) layer press(Pa)
  lyrprs_lyr*_0000_analinc - layer pressure increments in layer space (Pa, old in db)
  salint_lyr*_0000_analinc - salinity analysis increments in layer space (PSU)
  salint_lyr*_0000_analfld - updated salinity forecast in layer space (PSU)
  seatmp_lyr*_0000_analinc - temperature analysis increments in layer space 
                             (potential T deg C)
  seatmp_lyr*_0000_analfld - updated temperature forecast in layer space 
                             (potential T deg C)

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
#import netCDF4
import importlib
#import struct
#import pickle
#from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap
#from mpl_toolkits.basemap import Basemap, shiftgrid

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
import mod_read_ncoda as rncoda
#importlib.reload(rncoda)
from mod_read_hycom import read_grid_topo, read_hycom, read_topo
from mod_read_hycom import zz_zm_fromDP
from mod_utils_fig import bottom_text


rdate   = '2022061600'   # NCODA time 0 
#rdate  = '2022111700'
#pthbin = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_para7b/hycom/ncoda_archv_inc/'
pthbin = '/scratch2/NCEPDEV/marine/Jim.Cummings/zulema/'
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'

btx = 'plot_ncoda_lyrinc_input.py'

# get bacground date:
rbgrd = rncoda.anls2bckground_date(rdate)

IDM = 4500
JDM = 3298
KDM = 41
# Use either original file names dumped from NCODA or
# soft-links names created in ncoda_archv shell
# but check directories
#flnm = 'lyrprs_pre_1o{0}x{1}_{2}_0000_analinc'.format(IDM,JDM,rdate)

FLDS = {
"prsdat" : ["lyrprs", "fcstfld"],
"prsinc" : ["lyrprs", "analinc"],
"saldat" : ["salint", "analfld"],
"salinc" : ["salint", "analinc"],
"tmpdat" : ["seatmp", "analfld"],
"tmpinc" : ["seatmp", "analinc"]
}


get_topo = True
ftopo = 'regional.depth'
fgrid = 'regional.grid'



huge = 1.e20
rg   = 9806.
dbpa = 1.0198*rg  # NCODA conversion db to Pa 
itst = 2412
jtst = 1798

if get_topo:
#  HH = read_topo(pthgrid,ftopo,nn,mm)
  LON, LAT, HH = read_grid_topo(pthgrid,ftopo,fgrid)
  get_topo = False  



# Read layer pressure from HYCOM background (f/cast)
fldrd = 'prsdat'
aa    = FLDS.get(fldrd)
fld   = aa[0]
fend  = aa[1]
sznm  = '{0}x{1}'.format(IDM,JDM)
fhr   = 0
if fldrd == 'prsdat':
  fhr = 24
flnm = '{0}_lyr_1o{1}_{2}_00{3:02d}_{4}'.format(fld,sznm,rbgrd,fhr,fend)
#flnm = fld+'_lyr_1o'+sznm+'_'+rbgrd+'_0000_'+fend
fina = pthbin+flnm

print('Reading '+fina)
ZM = rncoda.read_ncoda_inc(fina,IDM,JDM,KDM)
ZM = -ZM/rg  

# Construct interface depths from ZM:
print('Constructing interface depths ZI ...')
ZI = np.zeros((KDM+1,JDM,IDM))
for kk in range(KDM):
  print('  k = {0}'.format(kk+1))
  dz = np.abs(ZM[kk,:,:]-ZI[kk,:,:])
  dmm = ZM[kk,:,:]-dz
  dmm = np.where(dmm<HH,HH,dmm)  # HH <0, dmm < 0
  dmm = np.where(dmm>0.,0,dmm)   # land = 0 m
  ZI[kk+1,:,:] = dmm

# ==========================
# Plot selected W-E profiles  
# selected fields
#
import plot_sect as psct
importlib.reload(psct)

clrmp = psct.clrmp_BlGrWhYlOrRd()
clrmp.set_bad(color = [0.7,0.7,0.7]) 

#FPLT = ["prsinc","saldat","salinc","tmpdat","tmpinc"]
FPLT = ["prsinc"]
#FPLT = ["salinc"]
#FPLT = ["tmpinc"]
#FPLT = ["saldat"]
#FPLT = ["tmpdat"]
nplt = len(FPLT)

#SCTP = ["GoM1"]
SCTP = ["SeaJpn"]
nsct = len(SCTP)
nfg  = 0
for ii in range(nplt):
  fldnm = FPLT[ii]
  aa = FLDS.get(fldnm)
  fld = aa[0]
  fend = aa[1]
  flnm = fld+'_lyr_1o'+sznm+'_'+rdate+'_0000_'+fend
  fina = pthbin+flnm

  A3d = np.array([])
  A3d = rncoda.read_ncoda_inc(fina,IDM,JDM,KDM)  

  print_pnt = True    
  if fldnm == "prsinc":
#    A3d = A3d*1.0198  # conver db to m
    A3d = A3d/rg
    rmin = -60.
    rmax = 60.
  elif fldnm == "salinc":
    rmin = -0.1
    rmax = 0.1
  elif fldnm == "tmpinc":
    rmin = -1.5
    rmax = 1.5
  elif fldnm == "saldat":
    clrmp = copy(plt.cm.rainbow) 
    clrmp.set_bad(color = [0.7,0.7,0.7])
    rmin = 34.5
    rmax = 36.5
  elif fldnm == "tmpdat":
    clrmp = copy(plt.cm.rainbow) 
    clrmp.set_bad(color = [0.7,0.7,0.7])
    rmin = 3.
    rmax = 28.
 
#    clrmp = copy(plt.cm.seismic)



  for ixx in range(nsct):
    nfg += 1
    xsct = SCTP[ixx]
    SCT = psct.Jsections()

    if xsct == "GoM1":
      jp0 = 72  # GoM corresponds the itest, jtest in ncoda_archv_inc
    elif xsct == "SeaJpn":
      jp0 = 55

    stl = ('{0} NCODA layerspace input {1}, {2}/{3}/{4}'.\
           format(xsct,flnm,rdate[0:4],rdate[4:6],rdate[6:8]))
    Hb, Xsct, Zsct, Asct = psct.plot_A3d_Jsct(xsct, A3d, ZI, HH, LON, LAT,\
                                fgnmb=nfg, stl=stl, sct_show=True, \
                                rmin=rmin, rmax=rmax, btx=btx, \
                                clrmp=clrmp, zstart=-1000., dcntr=1, jlrs=jp0)

    if print_pnt:
      psct.print_2D(Zsct[:,jp0],Asct[:,jp0],kend=KDM)




