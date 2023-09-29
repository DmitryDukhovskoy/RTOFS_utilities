"""
  Find Gulfstream north wall using definition of
  Halkin & Rossby, 1985 The Structure and Transport of the Gulf Stream at 73°W
  intersection of the 12C and 400m isobath

  NAVO Guld Stream north wall data:
  From Gwen: 
   On Hera, you can find the NAVO Gulf Stream north wall data files at /scratch2/NCEPDEV/ovp/Lichuan.Chen/DCOM/$YYYYMMDD/wtxtbul.  We use the gs*.sub to generate the plots at https://polar.ncep.noaa.gov/global/fronts/.  

   You can also find the data files on WCOSS2.  They are located at /lfs/h1/ops/prod/dcom/$YYYYMMDD/wtxtbul (e.g., /lfs/h1/ops/prod/dcom/20230407/wtxtbul/gs091nw.sub). 

# hour and output to plot wrt to rdate - f/cast date:
# -24, ..., 0, ..., 12
# for hr <=0, output from nowcast/hindcast "n-24" - initial after incup,
#    n-12, ..., n00 - initial for f/cast
# for hr >0: f01, ..., f12, ... f196 (if run that long) - forecasts


  Dmitry Dukhovskoy, NOAA/NCEP/EMC April-May 2023

"""
import os
import numpy as np
import sys
import importlib
import matplotlib
import matplotlib.pyplot as plt
import datetime

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hausdorff')

import mod_read_hycom as mhycom
import mod_misc1 as mmisc
import mod_hausdorff_distance as mmhd
import mod_utils_fig as mufig
import mod_utils as mutil
import mod_rtofs as mrtofs
importlib.reload(mmisc)
import mod_gulfstream as mgulf
importlib.reload(mgulf)
import mod_time as mtime
importlib.reload(mtime)

# Interpolation depth
z0 = -400.

expt    = 'paraD'
rdate0  = '20230419'
sfx     = 'n-24'

# Figure output directory:
f_figsave = False
f_intract = True  # False - no figures shown
pthfig  = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/' + expt + '/fig/'
pthscr  = '/scratch1/NCEPDEV/stmp2/Dmitry.Dukhovskoy/'
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'

if not f_intract:
  print('Interactive mode is OFF')
  matplotlib.use('Agg')
  plt.close('all')
  plt.ioff()
else:
  plt.ion()


pthhcm = pthscr + 'rtofs_{0}/run_diagn/rtofs.{1}/'.format(expt,rdate0)
fhcm   = 'rtofs_glo.t00z.' + sfx + '.archv'
fina   = pthhcm + fhcm + '.a'
finb   = pthhcm + fhcm + '.b'

if len(rdate0) > 0:
  YR     = int(rdate0[0:4])
  MM     = int(rdate0[4:6])
  DD     = int(rdate0[6:8])
#  yrday  = mmisc.date_yearday(YR,MM,DD)
  yrday  = mtime.rdate2jday(rdate0)
  dnmb0  = mtime.rdate2datenum(rdate0)

#
# Date of plotted fields:
if sfx == 'n-24':
  dnmbP = dnmb0-1
  hr = 0
elif sfx[0] == 'f':
  hr = int(sfx[1:])
  dnmbP = dnmb0+float(hr)/24.


dvP = mtime.datevec(dnmbP)
YRp = dvP[0]
MMp = dvP[1]
DDp = dvP[2]
HRp = dvP[3]


huge = 1.e20
rg   = 9806.
IDM, JDM, KDM = mhycom.hycom_dim(fina,finb)

ftopo = 'regional.depth'
fgrid = 'regional.grid'
LON, LAT, HH = mhycom.read_grid_topo(pthgrid,ftopo,fgrid)

# Read layer pressures:
dH = np.zeros((KDM,JDM,IDM))
for kk in range(1,KDM+1):
  F,nn,mm,lmm = mhycom.read_hycom(fina,finb,'thknss',rLayer=kk)

  F = F/rg
  F[np.where(F>huge)] = np.nan
  F[np.where(F<0.001)] = 0.
  dH[kk-1,:,:] = F

ZZ, ZM = mhycom.zz_zm_fromDP(dH)
kdm = ZM.shape[0]
jdm = ZM.shape[1]
idm = ZM.shape[2]

# Read T:
A3d  = np.array([])
fld  = 'temp'
for lvl in range (1,kdm+1):
  F,n1,m1,l1 = mhycom.read_hycom(fina,finb,fld,rLayer=lvl)
  F[np.where(F>huge)] = np.nan
  if A3d.size == 0:
    A3d = F.copy()
    A3d = np.expand_dims(A3d, axis=0)
  else:
    F = np.expand_dims(F, axis=0)
    A3d = np.append(A3d, F, axis=0)

# 2D Arrays with S/T above and below z0 for interpolation
print('Searching top/btm values for {0}'.format(z0))
Atop = np.zeros((jdm,idm))*np.nan
Abtm = np.zeros((jdm,idm))*np.nan
Ztop = np.zeros((jdm,idm))*np.nan
Zbtm = np.zeros((jdm,idm))*np.nan
for kk in range(kdm-1):
  zm1 = np.squeeze(ZM[kk,:,:])
  zm2 = np.squeeze(ZM[kk+1,:,:])
  if np.nanmin(zm2) > z0:
    continue
  if np.nanmax(zm1) < z0:
    continue

  print('lvl={0}  max(z1)/min(z2) = {1:5.1f}/{2:5.1f}'.\
         format(kk+1,np.nanmax(zm1),np.nanmin(zm2)))

  [J,I] = np.where((zm1 >= z0) & (zm2 < z0))
  Atop[J,I] = A3d[kk,J,I]
  Abtm[J,I] = A3d[kk+1,J,I]
  Ztop[J,I] = ZM[kk,J,I]
  Zbtm[J,I] = ZM[kk+1,J,I]

# linear interp
# Lagrange polynomial:
Aint = Atop*((z0-Zbtm)/(Ztop-Zbtm)) + \
       Abtm*((z0-Ztop)/(Zbtm-Ztop))


# Function to print mouse click event coordinates
def onclick(event):
   print([event.xdata, event.ydata])


# Define region of interest:
# For the HYCOM  0.08 Global grid !!! 
II = [2565, 2565, 2950, 2950]
JJ = [1809, 2190, 2190, 1809]

# Region to exclude shallow regions
ISH = [2562, 2569, 2574, 2579, 2581, 2609, 2721, 3012, 3012, 2562]
JSH = [1808, 1812, 1816, 1823, 1866, 1911, 1991, 2012, 1772, 1770]

X, Y = np.meshgrid(np.arange(IDM), np.arange(JDM))
MS, IM, JM    = mmisc.inpolygon_v2(X,Y,II,JJ)  # Gulfstream region
Msh, Ish, Jsh = mmisc.inpolygon_v2(X,Y,ISH,JSH)  # shelves to exclude

T12 = Aint.copy()
# For continuous contour 
# along the shelf fill shelf with cold water:
T12[np.where( (HH<0) & (HH>=z0) & (LAT>24.2) )] = -2.
T12[np.where(MS==0)] = np.nan
T12[np.where((Msh==1) & (HH>z0))] = np.nan
# Patches near Bahamas
T12[1808:1815, 2571:2586] = np.nan
TCNT = mgulf.derive_contour(T12, tz0=12.)


# Read NAVO Gulf Stream northwall coordinates:
# adjust # of days around rdate to look for NAVO
# fronts using missing=... 
importlib.reload(mrtofs)
importlib.reload(mgulf)
pthnavo = '/scratch2/NCEPDEV/ovp/Lichuan.Chen/DCOM/'
XNW_navo, YNW_navo, navonmb = mgulf.read_navogs(rdate0, pthnavo, missing=1)

print('Mapping NAVO lon/lat --> RTOFS index space ...')
INW_navo, JNW_navo = [], []
knavo = len(XNW_navo)
for ik in range(knavo):
  if ik%50 == 0:
    print('  {0:5.2f}% done ...'.format(float(ik)/float(knavo)*100.))

  xn0 = XNW_navo[ik]
  yn0 = YNW_navo[ik]

#  inavo, jnavo = mutil.interp_indx_lonlat(xn0, yn0, LON, LAT)
  inavo, jnavo = mrtofs.interp_lonlat2indx(xn0, yn0, LON, LAT)
  INW_navo.append(inavo)
  JNW_navo.append(jnavo)

# Truncate NAVO contour to match RTOFS if NAVO is outside the region
INW_navo = np.array(INW_navo)
JNW_navo = np.array(JNW_navo)
XNW_navo = np.array(XNW_navo)
YNW_navo = np.array(YNW_navo)

II   = np.array(II)
JJ   = np.array(JJ)
Pin  = mmisc.inpolygon(INW_navo,JNW_navo,II,JJ)
#Iout = np.where(Pin == False)[0]
Iin  = np.where(Pin == True)[0]
if len(Iin) < len(INW_navo):
  INW_navo = INW_navo[Iin]
  JNW_navo = JNW_navo[Iin]
  XNW_navo = XNW_navo[Iin]
  YNW_navo = YNW_navo[Iin]

#
# Compute Modified Hausdorff Distance
# Truncate RTOFS contour to NAVO
# Make sure RTOFS contour is longer than NAVO
# Next, map RTOFS Gulf Stream indices --> lon/lat space
# Compute MHD
INW_rtofs, JNW_rtofs = mgulf.adjust_rtofs_gsnw(TCNT[:,0], TCNT[:,1],\
                       INW_navo, JNW_navo, nskip=-1)

print('Mapping RTOFS GS northwall index ---> Lon/Lat space ...')
XNW_rtofs, YNW_rtofs = [], []
krt = len(INW_rtofs)
for ik in range(krt):
  if ik%50 == 0:
    print('  {0:5.2f}% done ...'.format(float(ik)/float(krt)*100.))

  xrt, yrt = mrtofs.interp_indx2lonlat(INW_rtofs[ik], JNW_rtofs[ik], LON, LAT)
  XNW_rtofs.append(xrt)
  YNW_rtofs.append(yrt)
              

XNW_rtofs = np.array(XNW_rtofs)
YNW_rtofs = np.array(YNW_rtofs)

knavo = len(XNW_navo)
P = np.zeros((krt,2))
Q = np.zeros((knavo,2))
P[:,0] = XNW_rtofs
P[:,1] = YNW_rtofs
Q[:,0] = XNW_navo
Q[:,1] = YNW_navo
mhdGS = mmhd.modifHD(P, Q, geo2cart=True)

print('=======   Start Plotting   =========')
from matplotlib import cm
#clrs   = cm.get_cmap('viridis',200)
clrs   = cm.get_cmap('rainbow',200)
clrs.set_bad(color=[0.3,0.3,0.3])
#  clrs.set_under(color=[0.8,0.8,0.8])
#  clrs.set_under(color=[0.8,0.7,1])

plt.ion()

rmin = 0.
rmax = 28.
tz0  = 12.0

Xt = TCNT[:,0]
Yt = TCNT[:,1]

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()

ax1 = plt.axes([0.1, 0.2, 0.7, 0.7])

Tplt = np.squeeze(A3d[0,:,:])
im1 = ax1.pcolormesh(Tplt, shading='flat', \
               cmap=clrs,\
               vmin=rmin, \
               vmax=rmax)
#  im1.set_clim(rmin,rmax)
#  ax1.contour(HH, [0.0], colors=[(0,0,0)], linewidths=1)
#  ax1.contour(T12, [tz0], colors=[(0.,0.4,0.6)], linewidths=1)
clr_rtofs = [0,0.8,1]
clr_navo  = [0,0,0]

# Plot Gulf Stream front:
ln1, = ax1.plot(Xt,Yt,'-',color=clr_rtofs, label="RTOFS")
# NAVO Gulf Stream North Wall:
ln2, = ax1.plot(INW_navo, JNW_navo,'-',color=clr_navo, label="NAVO")

lons = np.linspace(-180,180,73)
lats = np.linspace(-90,90,37)
ax1.contour(LON, lons, linestyles='dotted', colors=[(0.8,0.8,0.8)])
ax1.contour(LAT, lats, linestyles='dotted', colors=[(0.8,0.8,0.8)])

#
xlim1 = 2520
xlim2 = 3020
ylim1 = 1760
ylim2 = 2250 
#  ax1.axis('equal')
ax1.axis('scaled')
ax1.set_xlim([xlim1,xlim2])
ax1.set_ylim([ylim1,ylim2])
ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax1.set_xticks([])
ax1.set_yticks([])

# Actual limits:
Ylim = ax1.get_ylim()
Xlim = ax1.get_xlim()
Yl1  = int(Ylim[0])
Yl2  = int(Ylim[1])
Xl1  = int(Xlim[0])
Xl2  = int(Xlim[1])
# Put lon/lats on axis:
lon1 = LON[Yl1,Xl1]
lon2 = LON[Yl2,Xl2]
lat1 = LAT[Yl1,Xl1]
lat2 = LAT[Yl2,Xl2]

iLN1 = np.min(np.where(lons>=lon1)[0])
iLN2 = np.max(np.where(lons<=lon2)[0])
iLT1 = np.min(np.where(lats>=lat1)[0])
iLT2 = np.max(np.where(lats<=lat2)[0])
dltx = 1
dlty = 1
if iLN2-iLN1 >= 8:
  dltx = 2
if iLT2-iLT1 >= 8:
  dlty = 2

# X-axis labels
for ikk in range(iLN1,iLN2+1,dltx):
  xx0 = lons[ikk]
  yy0 = lat1   # for Mercator part of the grid, lat = const along j=fixed
  ii0, jj0 = mrtofs.find_indx_lonlat(xx0, yy0, LON, LAT)
  jj0 = Yl1-20
  xstl = '{0:3.1f}W'.format(abs(xx0))
  ax1.text(ii0, jj0, xstl, 
           fontsize=12,
           horizontalalignment='center')

# Y-axis labels
for ikk in range(iLT1,iLT2+1,dlty):
  yy0 = lats[ikk]
  xx0 = lon1
  ii0, jj0 = mrtofs.find_indx_lonlat(xx0, yy0, LON, LAT)
  ii0 = Xl1-10
  ystl = '{0:3.1f}N'.format(abs(yy0))
  if jj0 > Yl2:
    continue
  ax1.text(ii0, jj0, ystl, 
           fontsize=12,
           verticalalignment='center',
           horizontalalignment='right')

MHDinfo = 'Modif Hausdorff Dist RTOFS/NAVO = {0:5.1f} km'.format(mhdGS)
ctl = '{7} SST GS front 12C/400m, {5} {1}/{2:02d}/{3:02d}:{4:02d}\n {6}'.\
       format(z0, YR, MM, DD, hr, fhcm, MHDinfo, expt)
ax1.set_title(ctl)

ax1.set_xlim(Xlim)
ax1.set_ylim(Ylim)

# Select coordinate of the region of interest:
f_setrgn = False
if f_setrgn:
# Bind the button_press_event with the onclick() method
  fig1.canvas.mpl_connect('button_press_event', onclick)

ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
             ax1.get_position().y0,0.02,
             ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, extend='both')
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.tick_params(direction='in', length=12)

# Legend:
DV = mtime.datevec(navonmb)
ssinf = 'RTOFS: {0}/{1:02d}/{2:02d}:{3:02d}\n'.format(YR,MM,DD,hr)
ssinf = ssinf + 'NAVO: {0}/{1:02d}/{2:02d}:{3:02d}\n'.\
         format(DV[0],DV[1],DV[2],DV[3])
ax5 = plt.axes([0.7, 0.03, 0.2, 0.13])
lgd = plt.legend(handles=[ln1,ln2], loc='upper left')
ax5.text(0,0.01,ssinf)
ax5.axis('off')


btx = 'gulfstream_wall.py'
mufig.bottom_text(btx, pos=[0.1, 0.03])

if f_figsave:
  if not os.path.exists(pthfig):
    print('Creating ' + pthfig)
    os.makedirs(pthfig)

  fgtype = 'map'
  Regn   = 'NAtl'
  fld1   = 'gsfront'
#    fgnm   = 'GSnwall_' + expt + rdate + sfx + '.png' 
  fgnm   = '{0}_{1}_{2}_{3}_{4}_{5}.png'.format(fgtype,expt,rdate0,sfx,Regn,fld1)
  fpigout = pthfig + fgnm 
  print('Saving figure ---> ' + fpigout)
  plt.savefig(fpigout)

  


