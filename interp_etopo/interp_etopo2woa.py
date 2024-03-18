"""
  Interpolate 1/12.5 degree GLB topo to WOA grid
  Bilinear interpolation

  Dmitry Dukhovskoy, NOAA NESDIS NCEI
  June 2022
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import netCDF4
import importlib
from netCDF4 import Dataset as ncFile

sys.path.append('/data/home004/dmitry.dukhovskoy/python_codes/MyPython')
import mod_utils_fig
from mod_utils_fig import bottom_text
from mod_utils_fig import plot_fld2D

grd=0.25
if grd==0.25:
  cgrd=4
woa='woa18'
seas=13    # season: 1-12 - monthly, 13-16 - spring, summer, ...

pthout = '/data/ncei2/work/dmitry.dukhovskoy/topo/woa/'
pthhycom = '/data/ncei2/work/dmitry.dukhovskoy/topo/hycom_topo/'

urlT='https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/temperature/decav/{0}/'.\
      format(grd)
tfnm='{0}_decav_t{1:02d}_{2:02d}.nc'.format(woa,seas,cgrd)


def read_field(furl,varnm):
  print("Reading {1} from {0}".format(furl,varnm))
  nc=ncFile(furl)
# lookup a variable
  dmm0 = nc.variables[varnm][:].squeeze()
  dmm = np.copy(dmm0)
  return dmm

iz=0
tvar='t_an'

furl=urlT+tfnm
ZZ=read_field(furl,'depth')
ZZ=-abs(ZZ)
latW=read_field(furl,'lat')
lonW=read_field(furl,'lon')
jdmW=latW.shape[0]
idmW=lonW.shape[0]

#
# Make sure lon is -180,180
lonW=np.where(lonW>180.,lonW-360.,lonW)


#
# Read HYCOM topo and grid:
f_hycom=0
if f_hycom==1:
  ftopo = pthhycom+'depth_GLBa0.08_09.nc'
#nc = ncFile(ftopo)
  LATH = read_field(ftopo,'Latitude')
  LONH = read_field(ftopo,'Longitude')
  HH   = read_field(ftopo,'bathymetry')
  HH[np.where(HH>1.e20)] = np.nan
  HH = -HH
  HH[np.where(np.isnan(HH))]=100.  # land
  jdm  = HH.shape[0]
  idm  = HH.shape[1]
  Lmin = np.min(LATH)

fland = 100.
HWOA = np.zeros((jdmW,idmW))+1.e20

Iall = np.where(HWOA.flatten()>=0)[0]  
nall = Iall.shape[0]

from mod_misc1 import dist_sphcrd


#
# Read binary ETOPO2
pthetopo = '/data/ncei3/nodc/data/oc5.ocl/masks/'
#fbin = 'etopo2real.d'
fbin = 'etopo2.real'
fetopo = pthetopo+fbin
fnc = pthetopo+'etopo2.grd'

#LATE = read_field(fnc,'y_range')
#LONE = read_field(fnc,'x_range')
#HE = read_field(fnc,'z')

# Read binary
#dt = np.dtype('f8')
f_bin=0
if f_bin==1:
  idmE = 10800
  jdmE = 5400
  fid = open(fetopo,'rb')
  data = np.fromfile(fid,dtype='<i4', count=idmE*jdmE)
#data = np.fromfile(fid,dtype='<f4', count=idmE*jdmE)
  fid.close()

  HE = data.reshape((jdmE,idmE))
  HE[np.where(HE==-1)]=100.

lonE = np.arange(-180.,180.,1./30.,dtype=float)
latE = np.arange(-90,90,1./30.,dtype=float)

idmE = 10801
jdmE = 5401
data = np.zeros((idmE*jdmE))
lon  = np.zeros((idmE))
lat  = np.zeros((jdmE))
f_ascii=1
if f_ascii == 1:
  print('Reading ascii etopo2')
  HE = np.zeros((jdmE,idmE))
# Invert lat direction: 1=S. Pole, Last = N. Pole
# And lon is flipped E-W - need to reorient
  fetopo = pthetopo+'etopo2.ascii'
  fid = open(fetopo,'r')
#  data = np.genfromtxt(fid, delimiter=' ')
  cc=0
  icc=idmE
  jcc=0
  for line in fid:
    aa=line.strip()
    cll = line.split() 
    dmm = float(cll[2])
    x0  = float(cll[0])
    y0  = float(cll[1])

#    data[cc]=dmm
    cc += 1
    icc += 1
    if icc >= idmE:
      icc = 0
      jcc += 1

    if jcc == 1:
      lon[icc]=x0

    if icc == 0:
      lat[-jcc]=y0
 
    HE[-jcc,icc] = dmm

    if cc % 1e6 == 0:
      print('cc={0}: lon={1} lat={2} depth={3}  {4:3.1f}% done'.\
            format(cc,x0,y0,dmm,cc/(idmE*jdmE)*100))

  fid.close()
  
  lonE = lon
  latE = lat   


# Make sure lon is -180,180
lonE=np.where(lonE>180.,lonE-360.,lonE)

#nbytes = idmE*jdmE*4
#with open(fetopo,"rb") as ff:
#  ff.seek(0,0)
#  f1 = ff.read(1)

#ff=open(fetopo,"rb")
#aa=ff.read(1)
#ff.close()


from copy import copy
import matplotlib.colors as colors
f_chck=0
if f_chck==1:
  clrmp = copy(plt.cm.jet)
  clrmp.set_over(color=[0.5, 0.5, 0.5])


  plt.ion()
  fig1 = plt.figure(1, figsize=(8,8), constrained_layout=False)
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  cl1=-8000.
  cl2=0.
  f_indx=1
  if f_indx>0:
    im1 = ax1.pcolormesh(HE,shading='flat',cmap=clrmp,
                norm=colors.Normalize(vmin=cl1, vmax=cl2))
    ax1.axis('scaled')
    ax1.set(xlim=(0,idmE),ylim=(0,jdmE))
  else:
    im1 = plt.pcolormesh(lonE,latE,HE,shading='nearest',cmap=clrmp,
                 norm=colors.Normalize(vmin=cl1, vmax=cl2))
    ax1.axis('scaled')
    ax1.set(xlim=(-180,180),ylim=(-90,90))

  ax1.autoscale(enable=True, axis='both', tight=True)

  fig1.colorbar(im1,ax=ax1,orientation='horizontal')

  im1.set_clim(cl1,cl2)

#  plt.colorbar(im1)
  plt.show() 


#
import mod_bilinear
#importlib.reload(mod_bilinear)
from mod_bilinear import basisFn_RectRef
from mod_bilinear import phi_x0y0
from mod_bilinear import map_x2xhat
from mod_bilinear import phi_x0y0
from mod_bilinear import bilin_interp
#
# Find basis functions for a reference rectangle:
phi1,phi2,phi3,phi4 = basisFn_RectRef()


import timeit
kcc = 0
tic = timeit.default_timer()
ticR = timeit.default_timer()
print('Interpolating HYCOM topo -----> WOA grid')

# Check that min/max lon same sign
#if np.min(lonE):
#for ikk in range(1):  # testing
for ikk in range(nall):
  I1 = Iall[ikk]
  jj, ii = np.unravel_index(I1,HWOA.shape)
  kcc += 1

  if (kcc % 20000) == 0:
    toc = timeit.default_timer()
    try:
      Hmin = np.nanmin(HWOA[np.where(HWOA<=0.)])
      Hmax = np.nanmax(HWOA[np.where(HWOA<=0.)])
    except:
      print(' Hmin/Hmax are all land ')
      Hmin = np.nanmin(HWOA[np.where(HWOA<=1.e10)])
      Hmax = np.nanmax(HWOA[np.where(HWOA<=1.e10)])


    print(' {0:5.2f}% done {1:6.2f} min tot, {2:6.2f} min, lat={3:5.1f} ...'.\
            format(kcc/nall*100,(toc-tic)/60,(toc-ticR)/60,lat[jj]))
    print('  Min/Max Topo WOA = {0:6.2f} / {1:6.2f} m'.format(Hmin,Hmax))
    print(' ')

    ticR = timeit.default_timer()

  x0 = lonW[ii]
  y0 = latW[jj]
#  if y0 < Lmin:
#    continue
  ii1=np.where(lonE<=x0)[0][-1]
  ii2=np.where(lonE>=x0)[0][0]
  jj1=np.where(latE<=y0)[0][-1]
  jj2=np.where(latE>y0)[0][0]

  if ii2-ii1 != 1:
    print('*** ERR: check indices ii1={0}, ii2={1}'.format(ii1,ii2))
    sys.exit()
  
  if jj2-jj1 != 1:
    print('*** ERR: check indices ii1={0}, ii2={1}'.format(jj1,jj2))
    sys.exit()

# Check boundary:
  if ii1 == ii2:  # case when pnt on the 180W
    ii1 = ii2-1
  
  if jj1 == jj2:   # pnt on 90N
    jj1 = jj2-1

  XX = np.array([lonE[ii1],lonE[ii2],lonE[ii2],lonE[ii1]], dtype=float)
  YY = np.array([latE[jj1],latE[jj1],latE[jj2],latE[jj2]], dtype=float)
  HT = np.array([HE[jj1,ii1], \
                 HE[jj1,ii2], \
                 HE[jj2,ii2], \
                 HE[jj2,ii1]])

#  DS = dist_sphcrd(y0,x0,LATH,LONH)
#
# Map X,Y ---> Xhat, Yhat on reference quadrialteral 
# i.e. map WOA grid coordinate to a reference quadrilateral 
# to do bilinear interpolation 
  xht, yht = map_x2xhat(XX,YY,x0,y0)

# Perform interpolation on reference rectangle, that is 
# similar to interp on actual rectangle
  hintp = bilin_interp(phi1,phi2,phi3,phi4,xht,yht,HT)

#  if hintp > 0.:
#    HWOA[jj,ii] = fland
#  else:
#    HWOA[jj,ii] = hintp
  HWOA[jj,ii] = hintp


print('Interplation finished ')

# Saving output:
fout = pthout+'etopo2woa025_bilinear.dat'
print(' Saving interpolated TOPO to WOA 0.25 grid: {0}'.format(fout))
pf1 = open(fout,'wb')
HWOA.astype('f').tofile(pf1)
pf1.close()

# Lon/lat:
fout = pthout+'lonWOA025.dat'
print(' Saving Lon WOA 0.25 grid: {0}'.format(fout))
pf1 = open(fout,'wb')
lonW.astype('f').tofile(pf1)
pf1.close()

fout = pthout+'latWOA025.dat'
print(' Saving Lat WOA 0.25 grid: {0}'.format(fout))
pf1 = open(fout,'wb')
latW.astype('f').tofile(pf1)
pf1.close()



print(' =========  Done  ========')
print(' ')


btx = 'interp_hycom2woa.py'

f_chck=0
if f_chck==1:
  clrmp = copy(plt.cm.jet)
  clrmp.set_over(color=[0.5, 0.5, 0.5])


  plt.ion()
  fig1 = plt.figure(2, figsize=(8,8), constrained_layout=False)
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  cl1=-8000.
  cl2=0.
  f_indx=1
  if f_indx>0:
    im1 = ax1.pcolormesh(HWOA,shading='flat',cmap=clrmp,
                norm=colors.Normalize(vmin=cl1, vmax=cl2))
    ax1.axis('scaled')
    ax1.set(xlim=(0,idmW),ylim=(0,jdmW))
  else:
    im1 = plt.pcolormesh(lonW,latW,HWOA,shading='nearest',cmap=clrmp,
                 norm=colors.Normalize(vmin=cl1, vmax=cl2))
    ax1.axis('scaled')
    ax1.set(xlim=(-180,180),ylim=(-90,90))

  ax1.set_title('TOPO WOA0.25 from ETOPO2, bilinear interp')
  ax1.autoscale(enable=True, axis='both', tight=True)

  fig1.colorbar(im1,ax=ax1,orientation='horizontal')

  im1.set_clim(cl1,cl2)

  bottom_text(btx,pos=[0.1,0.1])
#  plt.colorbar(im1)
  plt.show()







