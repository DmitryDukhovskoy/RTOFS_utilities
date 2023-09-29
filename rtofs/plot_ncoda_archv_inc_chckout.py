"""
  Rho background and rho from (T+Tincr,S+Sincr)
  Plot diagnostic output fields from
  ncoda_archv_inc.f / ncoda_arch_lyrinc.f code
  to check increments/pressure/density changes
  caused by ncoda incr etc
  These are intermediate fields - not final for checking only
  Note that lyr pressure change is the one caused by 
  rho(T+incrT, S+incrS) and not necessarily the actual
  lr pr change that will be saved at the end
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import netCDF4
import importlib
import struct
from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap
#from mpl_toolkits.basemap import Basemap, shiftgrid
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
import mod_read_ncoda as rncoda

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')

#rdate   = '20221117'
rdate   = '20220616'
pth1    = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_para7b/hycom/ncoda_archv_inc'

pthbin  = pth1+'/'
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'

ncoda_levels = False          # old version of ncoda_archv_inc saves on NCODA z-levels
                              # updated version - on HYCOM layers 
                              # these intermediate fields are at mid-layer points
frhobg    = 'rho_bckgr'       # Rho of background field on NCODA z-levels 
frhoincr  = 'rho_added_incr'  # Rho (T+Tincr,S+Sincr)
fprhoincr = 'p_rho_incr'      # Pressure change (lr depth)

IDM = 4500
JDM = 3298
KDM = 41

# ==========
# pltrho   - Plot background rho on NCODA z-levels
# pltrhoi  - Plot rho computed from T,S after adding increments
# pltdrho  - Plot difference between the two rho fields
# pltdplr  - Plot pressure increments (layer interface displacement)
#            from density change after adding increments
f_pltrho  = True
f_pltrhoi = False
f_pltdrho = False
f_pltdplr = False


# Plot vertical section of background and increm rho
# and difference btw the two
SCTP = ["GoM1","SeaJpn"]
nsct = len(SCTP)


get_topo = True
ftopo = 'regional.depth'
fgrid = 'regional.grid'

import mod_read_hycom
#importlib.reload(mod_read_hycom)
from mod_read_hycom import read_grid_topo, read_hycom, \
                           read_topo, read_hycom_restart
import mod_utils as utls
#importlib.reload(utls)
from mod_read_hycom import read_2Dbinary
from mod_read_hycom import zz_zm_fromDP
from mod_utils_fig import bottom_text

btx = 'plot_ncoda_archv_inc_chckout.py'

def print_1D(A,wd=8,prc=2):
  ndim1 = A.shape[0]
  for k in range (ndim1):
    print('{0}: {1:{width}.{precis}f}'.format(k+1,A[k],width=wd,precis=prc))

get_topo = True
if get_topo:
#  HH = read_topo(pthgrid,ftopo,nn,mm)
  LON, LAT, HH = read_grid_topo(pthgrid,ftopo,fgrid)
  get_topo = False  

huge = 1.e20
rg   = 9806.

# Background rho
fina = pthbin+frhobg+'.a'
finb = ""      # only binary
print('Processing '+fina)

Rho_bg = np.array([])
for kk in range (1,KDM+1):
  F = read_2Dbinary(fina,IDM,JDM,kk)
  F[np.where(F>huge)] = np.nan
  if Rho_bg.size == 0:
    Rho_bg = F.copy()
    Rho_bg = np.expand_dims(Rho_bg, axis=0)
  else:
    F = np.expand_dims(F, axis=0)
    Rho_bg = np.append(Rho_bg, F, axis=0)


# Rho from T+Tincr,S+Sincr
fina = pthbin+frhoincr+'.a'
finb = ""      # only binary
print('Processing '+fina)

Rho_incr = np.zeros((KDM,JDM,IDM))
for kk in range (1,KDM+1):
  F = read_2Dbinary(fina,IDM,JDM,kk)
  F[np.where(F>huge)] = np.nan
  Rho_incr[kk-1,:,:] = F


# Pressure change from rho difference
# layer interf depth change
fina = pthbin+fprhoincr+'.a'
find = ""
print('Processing '+fina)

dHlr = np.zeros((KDM,JDM,IDM))
for kk in range(1,KDM+1):
  F = read_2Dbinary(fina,IDM,JDM,kk)
  F[np.where(F>huge)] = np.nan
  dHlr[kk-1,:,:] = F/rg       # m



# Read NCODA interface & mid-pnt depths:
# get bacground date:
rbgrd = rncoda.anls2bckground_date(rdate)
if ncoda_levels:
  ZI, kzi, ZM, kzm = utls.ncoda_depths()
else:
# Read from b/ground archv or lyr prs ncoda input:
  aa    = "prsdat"
  fld   = "lyrprs"
  fend  = "fcstfld"
  sznm  = '{0}x{1}'.format(IDM,JDM)
  fhr   = 0
  fhr   = 24
  flnm  = '{0}_lyr_1o{1}_{2}_00{3:02d}_{4}'.format(fld,sznm,rbgrd,fhr,fend)
  #flnm = fld+'_lyr_1o'+sznm+'_'+rbgrd+'_0000_'+fend
  fina  = pthbin+flnm

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



# ================================
def find_indx_lonlat(x0,y0):
  """
  Find closest grid point to lon/lat coordinate
  """
  if x0 > 180.:
    x0 = x0-360.
  
  dmm = np.sqrt((LON-x0)**2+(LAT-y0)**2)
  jj0, ii0 = np.where(dmm == np.min(dmm))

  return ii0[0], jj0[0]

# ================================
kdm = Rho_incr.shape[0]
jdm = Rho_incr.shape[1]
idm = Rho_incr.shape[2]



# ==========================
# Plot selected W-E profiles  
import plot_sect as psct
importlib.reload(psct)

if f_pltrho:
  nfg = 0
  for ii in range(nsct):
    nfg += 1
    xsct = SCTP[ii]
    SCT = psct.Jsections()

    stl = ('ncoda_archv_inc: Rho_bgr {0} {1}/{2}/{3}'.\
           format(xsct,rdate[0:4],rdate[4:6],rdate[6:]))
    Hb, XX, Asct = psct.plot_rho_Jsct(xsct, Rho_bg, ZI, HH, LON, LAT,\
                                fgnmb=nfg, stl=stl, sct_show=True, \
                                rmin=32.2, rmax=37.2, btx=btx)


if f_pltrhoi:
  nfg = 10
  for ii in range(nsct):
    nfg += 1
    xsct = SCTP[ii]
    SCT = psct.Jsections()

    stl = ('ncoda_archv_inc: Rho_updated {0} {1}/{2}/{3}'.\
           format(xsct,rdate[0:4],rdate[4:6],rdate[6:]))
    Hb, XX, Asct = psct.plot_rho_Jsct(xsct, Rho_incr, ZI, HH, LON, LAT,\
                                fgnmb=nfg, stl=stl, sct_show=True, btx=btx)

if f_pltdrho:
  nfg = 20
  for ii in range(nsct):
    nfg += 1
    xsct = SCTP[ii]
    SCT = psct.Jsections()

    clrmp = copy(plt.cm.BrBG_r)
#    dRho = Rho_bg-Rho_incr
    dRho = Rho_incr-Rho_bg  
    stl = ('ncoda_archv_inc: (Rho_updt-Rho_bgr) {0} {1}/{2}/{3}'.\
           format(xsct,rdate[0:4],rdate[4:6],rdate[6:]))
    Hb, XX, Asct = psct.plot_rho_Jsct(xsct, dRho, ZI, HH, LON, LAT,\
                   fgnmb=nfg, stl=stl, sct_show=True, btx=btx, \
                   rmin=-0.25, rmax=0.25, clrmp=clrmp)

if f_pltdplr:
  nfg = 30
  for ii in range(nsct):
    nfg += 1
    xsct = SCTP[ii]
    SCT = psct.Jsections()

    clrmp = copy(plt.cm.RdBu_r)
    stl = ('ncoda_archv_inc: dP(m)  {0} {1}/{2}/{3}'.\
           format(xsct,rdate[0:4],rdate[4:6],rdate[6:]))
    Hb, XX, Asct = psct.plot_rho_Jsct(xsct, dHlr, ZI, HH, LON, LAT,\
                fgnmb=nfg, stl=stl, sct_show=True, btx=btx, \
                rmin=[], rmax=[], clrmp=clrmp)



# Plot 2D rho map for specified depth level
f_pltrdp = False
if f_pltrdp:
  fgnmb1 = 1
  print('Plotting rho ...')
  clrmp = copy(plt.cm.afmhot_r)
  clrmp.set_bad(color=[0.7,0.7,0.7])
  plt.ion()
  fig1 = plt.figure(fgnmb1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  im1 = plt.pcolormesh(F, vmin=0, vmax=40, cmap=clrmp)
  plt.contour(HH,[0.0],colors=[(0,0,0)],linewidths=1)
 
  ax2 = fig1.add_axes([ax1.get_position().x1+0.02, 
               ax1.get_position().y0,0.02,
               ax1.get_position().height])
  clb = plt.colorbar(im1, cax=ax2, extend='max')
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  clb.ax.set_yticklabels(ticklabs,fontsize=10)

  plt.sca(ax1)

  ctl = 'Inverse min(lr_thkns/dH) wrt {0:3.2f}, {1}'.\
        format(rdp0,rdate+'/'+flnm)
  ax1.set_title(ctl)
 
  bottom_text(btx)




