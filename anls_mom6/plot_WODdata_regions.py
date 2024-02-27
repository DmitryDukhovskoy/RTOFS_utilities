"""
  Plot T/S Profile Observations - downloaded from WOD18 website
  https://www.ncei.noaa.gov/access/world-ocean-database/bin/getwodyearlydata.pl
  show on orthonormal map

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
import datetime
import pickle
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from netCDF4 import Dataset as ncFile

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

from mod_utils_fig import bottom_text
import mod_mom6_valid as mom6vld
import mod_read_hycom as mhycom
import mod_misc1 as mmisc
import mod_time as mtime
import mod_WODdata as mwod
import mod_mom6 as mom6util
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.colors as colors
import matplotlib.mlab as mlab
importlib.reload(mwod)


# WOD data
pthwod = '/scratch1/NCEPDEV/stmp4/Dmitry.Dukhovskoy/WOD_profiles/'
pthoutp= '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data_anls/MOM6_CICE6/ts_prof/'
pthrun = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/MOM6_run/' + \
         '008mom6cice6_003/'


# Select lon, lat to search for WOD profiles
YR1    = 2018
YR2    = 2022
mo1    = 1
mo2    = 12
#regn   = 'AmundsAO'
#regn   = 'NansenAO'
#regn   = 'MakarovAO'
#regn  = 'CanadaAO'

REGNS = mom6vld.ts_prof_regions()
btx   = 'plot_WODdata_regions.py'

clr_drb = [0.8, 0.4, 0]
clr_pfl = [0, 0.4, 0.8]
clr_ctd = [1, 0, 0.8]

clr_amd = [0., 0.3, 0.52]
clr_nns = [0.75, 0.25, 0.25]
clr_mkr = [0.45, .72,  0.54]
clr_cnd = [1.,  0.85, 0.4]

lon0=-10
lat0=70


pthgrid   = pthrun + 'INPUT/'
fgrd_mom  = pthgrid + 'regional.mom6.nc'
ftopo_mom = pthgrid + 'ocean_topog.nc'
LON, LAT  = mom6util.read_mom6grid(fgrd_mom, grdpnt='hpnt')
HH        = mom6util.read_mom6depth(ftopo_mom)


f_new = True
for regn in REGNS:
  x0 = REGNS[regn]["lon0"]
  y0 = REGNS[regn]["lat0"] 
  dx = REGNS[regn]["dlon"]
  dy = REGNS[regn]["dlat"]

  print(f'\nPlotting WOD for {regn}\n')

  fuid_out = f'uidWOD_{regn}_{YR1}-{YR2}.pkl'
  dflout   = os.path.join(pthoutp,fuid_out)

# Load saved UID:
  print(f'Loading {dflout}')
  with open(dflout, 'rb') as fid:
    [SIDctd, SIDdrb, SIDpfl] = pickle.load(fid)

  Xctd, Yctd = mwod.derive_lonlatWOD(SIDctd, pthwod, 'CTD')
  Xdrb, Ydrb = mwod.derive_lonlatWOD(SIDdrb, pthwod, 'DRB')
  Xpfl, Ypfl = mwod.derive_lonlatWOD(SIDpfl, pthwod, 'PFL')


  plt.ion()
  obtype = 'DRB'

  if f_new:
    m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution='l')
    xh, yh = m(LON,LAT)  # modl grid coordinates on the projections
    fig1 = plt.figure(1,figsize=(9,9))
    plt.clf()
    ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
    m.drawcoastlines()
    m.fillcontinents(color=(0.2,0.2,0.2))
    m.drawparallels(np.arange(-90.,120.,10.))
    m.drawmeridians(np.arange(-180.,180.,10.))
    m.contour(xh, yh, HH, [-4000,-3000,-2000,-1000],
              colors=[(0.8,0.8,0.8)],
              linestyles='solid')
    f_new = False


  xpfl_m, ypfl_m = m(Xpfl,Ypfl)  
  xdrb_m, ydrb_m = m(Xdrb,Ydrb)  
  xctd_m, yctd_m = m(Xctd,Yctd)  
  m.plot(xpfl_m,ypfl_m,'.', color=clr_pfl)
  m.plot(xdrb_m,ydrb_m,'.', color=clr_drb)
  m.plot(xctd_m,yctd_m,'.', color=clr_ctd)

# Plot specified region:
  if regn[:4] == 'Cana':
    clrrg = clr_cnd
  elif regn[:4] == 'Nans':
    clrrg = clr_nns
  elif regn[:4] == 'Amun':
    clrrg = clr_amd
  elif regn[:4] == 'Maka':
    clrrg = clr_mkr  

  ddx = 5
  nstps = int(2*dx/ddx) + 1
  XR  = [x0-dx]
  YR  = [y0-dy]
  for icc in range(nstps):
    xnxt = (x0-dx) + icc*ddx
    if xnxt > x0+dx:
      break
    XR.append(xnxt)
    YR.append(y0+dy)

  for icc in range(nstps):
    xnxt = (x0+dx) - icc*ddx
    if xnxt < x0-dx:
      break
    XR.append(xnxt)
    YR.append(y0-dy)

  XRp, YRp = m(XR,YR)
  m.plot(XRp,YRp,'-', color=(tuple(clrrg)))

sttl = f'Regions and  WOD obs {YR1}-{YR2} {mo1}-{mo2}'
ax1.set_title(sttl)

ax2 = plt.axes([0.78,0.82,0.18,0.15])
ax2.plot(0, 0.9, '.', ms=14, color=clr_ctd)
ax2.text(0.07, 0.9, 'CTD')
ax2.plot(0, 0.5, '.', ms=14, color=clr_drb)
ax2.text(0.07, 0.5, 'DRB')
ax2.plot(0, 0.1, '.', ms=14, color=clr_pfl)
ax2.text(0.07, 0.1, 'PFL')
ax2.set_xlim([-0.1, 0.4])
ax2.set_ylim([0, 2])
ax2.axis('off')

bottom_text(btx)


