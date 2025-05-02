"""
  Create inverse relaxation e-folding time scale
  Time scale can vary spatially allowing different relaxation
  rates for different parts of the domain

  Relaxation field is created using complex transformation/mapping technique
 
"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
import matplotlib.pyplot as plt
from yaml import safe_load

PPTHN = '/home/Dmitry.Dukhovskoy/python'
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')
import mod_time as mtime
import mod_mom6 as mmom6
import mod_utils as mutil
import mod_colormaps as mclrmps
import mod_misc1 as mmisc
from mod_utils_fig import bottom_text

# Select max and min relaxation time scales, hrs
# max relaxation - strongest, typically along the OBs
# min relaxation - somewhere in the domain where sea ice presents
#                  e.g. Bering strait
# rx_min, ry_min - approximate # of i, j pnts from the ice OBs
#                  i.e., from i=imax to Ber. Str. (342-200)
# relaxation time scales will be going to 0 away from the ice OBs
rate_max_hrs  = 24.                     # max relaxation time, hrs
Irate_max_sec = 1./(rate_max_hrs*3600.)  # relaxation rate, s-1

f_save    = False         # Save netcdf relax file
check_rlx = True         # Plot relaxation field
check_ref_domain = True  # Plot transformations of the reference domain

rlx_name = 'relax_rate' # name of the variable, should be the same in the SIS_input

btx  = 'relax_timescale.py'

if not f_save:
  print(f'WARNING: relaxation field is not saved, f_save: {f_save}')

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

# MOM6 NEP topo/grid:
run_name   = 'seasonal_fcst_daily'
pthtopo    = gridfls['MOM6_NEP'][run_name]['pthgrid']
fgrid      = gridfls['MOM6_NEP'][run_name]['fgrid']
ftopo_mom  = gridfls["MOM6_NEP"][run_name]["ftopo"]
outdir     = gridfls['MOM6_NEP'][run_name]['pthoutp']
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
# Hgrid lon. lat:
hlon, hlat  = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

# Define domain in a complex plane with 
# North-East boundaries along the Im axis (strongest relaxation)
# and slowest relax. rate at a distance dist_rlx 
# rpolar  = approximate # of grid points from the NE corner  of the NEP domain
# to ~55N
# Define a rectangular domain D, transformed relax zone is a square = rpolar
rpolarN = 250   # this defines the size of the relaxation zone along N. bndry axis (index space)
rpolarE = 260   # -"- -"- along the E. bndry axis
#rpolar = 240
Y = np.arange(-rpolarN,rpolarE+1)
X = np.arange(0,idm)
XR = np.arange(-X[-1],1)     # for rotated domain that is in the X<0 plane
YR = np.arange(-Y[-1],-Y[0]+1)  # rotated domain 
idmD = len(X)
jdmD = len(Y)

# Specify the number of grid points from the boundary to keep max relaxation
# Define a reference rectangular domain with relaxation decaying from the left bndry to the right
# Decrease the width of the region of strong relaxation along the N. bndry going south
# and increase the width along the E. bndry going south to have >0 relax over Eastern Bering shelf
irmax = 90    # where to begin exponential decay of the relaxation
ieN = 135     # start pnt (in domain D indices) where the width of max rlx starts changing
isN = 10 
irmaxN = 1    # width of max rlx at the S. end of N. boundary
isE = jdmD-230
ieE = jdmD
irmaxE = 120

# Adjust exponential decay of the relaxation off the Y=0 axis:
f_adj = False
sgmx = idm/4  # controls exponential decay in the Gaussian, decrease the denom. to slow the decay
RLX = np.zeros((jdmD, idmD))
for jj in range(jdmD):
  irmax0 = irmax
  if f_adj:
    if jj>=isN and jj<=ieN:
      irmax0 = irmaxN + int((irmax-irmaxN)/(ieN-isN-1)*(jj-isN))
    elif jj<isN:
      irmax0 = irmaxN
    elif jj>=isE and jj<=ieE:
      irmax0 = irmax + int((irmaxE-irmax)/(ieE-isE-1)*(jj-isE))

  aa = np.exp(-((X-irmax0)**2/sgmx**2))
  aa = aa/np.max(aa)
  aa[:irmax0] = aa[irmax0]
  RLX[jj,:] = aa*Irate_max_sec

# Perform 1st mapping using z**rexp
# Use the fact that mapped domain is symmetric wrt real axis X
rdnm = 1.4
rexp = 1/rdnm
RMAP1 = np.zeros((jdmD, idmD))*np.nan
RMAP2 = RMAP1.copy()*np.nan
jD0 = np.argmin(np.abs(Y))
#theta_dgr = 225. # rotation angle, sign in math. sense
theta_dgr = 270.-90./rdnm # rotation angle, sign in math. sense
for ii in range(idmD):
  for jj in range(jdmD):
    # Y Distance wrt to jD0:
    xx = X[ii]
    yy = Y[jj]
    RR  = np.sqrt(xx**2+yy**2)
    phi = np.arctan2(yy,xx)
    # Note that under the 1st mapping, RR should be (RR)**(1/rndm)
    # To keep RR unchanged during the transformation and keep the same grid dim
    # apply another transformation z=z*|z|**((rdnm-1)/rdnm), i.e. (RR)**(1/rdnm) --> RR
    x_map = RR*np.cos(phi/rdnm)
    y_map = RR*np.sin(phi/rdnm)
    imap = np.max(np.where(X<=x_map)[0])
    jmap = np.max(np.where(Y<=y_map)[0])
    #RMAP1[jmap, imap] = RLX[jj,ii]      # this may leave empty cells due to round off errors
    imm1 = np.max([imap-1,0])
    imm2 = np.min([imap+1,idmD])
    jmm1 = np.max([jmap-1,0])
    jmm2 = np.min([jmap+1,jdmD])
    RMAP1[jmm1:jmm2, imm1:imm2] = RLX[jj,ii]

    # Perform 2nd mapping - rotation by theta- degree angle:
    theta = theta_dgr*np.pi/180.
    RRmap1 = np.sqrt((x_map)**2 + (y_map)**2)
    phi1 = np.arctan2(y_map,x_map)
    rmap0 = RMAP2[jD0,jD0]    # debugging
    x_rot = RRmap1*np.cos(phi1+theta)
    y_rot = RRmap1*np.sin(phi1+theta)
    indx = np.where(XR>=x_rot)[0]
    jndx = np.where(YR>=y_rot)[0]
    if indx.size > 0 and jndx.size > 0:
      irot = np.min(np.where(XR>=x_rot)[0])
      jrot = np.min(np.where(YR>=y_rot)[0])
      #RMAP2[jrot,irot] = RLX[jj,ii]
      # To avoid gaps: fill neighboring grid cells
      ip1 = np.min([irot+1, idmD])
      im1 = np.max([irot-1, 0])
      jp1 = np.min([jrot+1, jdmD])
      jm1 = np.max([jrot-1, 0])
      RMAP2[jm1:jp1,im1:ip1] = RLX[jj,ii]

# Imbed transformed relax. field into relaxation array PSI
jdmR, idmR = RMAP2.shape
# The domain is in the 3rd quarter:
# domain: Y<=0 and x<=0
# subset part of the domain that contains transformed field
iyax0 = max(np.where(YR<=0)[0])
ixax0 = max(np.where(XR<=0)[0])
AA = RMAP2[:iyax0+1,:ixax0]    
AA = np.where(np.isnan(AA), 0., AA)
jdmA, idmA = AA.shape
isD = idm-idmA
jsD = jdm-jdmA
RLXIS = np.zeros((jdm,idm))
RLXIS[jsD:,isD:] = AA    # relaxation rate, s-1
# No relaxation in the Gulf of Alaska:
RLXIS[:575,170:] = 0.0


# Add land mask and southern domains = 0
lat_cut = 53.
RLXIS = np.where(HH>=0, 0.0, RLXIS)
RLXIS = np.where(hlat<lat_cut, 0.0, RLXIS)

# No relaxation in the G. Alaska:
RLXIS[:574,140:] = 0.0
RLXIS[:580,90:141] = 0.0

# For checking, relaxation time, hrs:
RLXHR = RLXIS.copy()
RLXHR = np.where(RLXHR==0., np.nan, RLXHR)
RLXHR = 1./RLXHR * 1/3600.

# Write relax time scale:
dflstat = os.path.join(pthtopo, 'ocean_static.nc')
ds_stat = xarray.open_dataset(dflstat)
ds_rlx  = ds_stat['wet'].copy()
ds_rlx.name = rlx_name
ds_rlx *= RLXIS
ds_rlx = ds_rlx.to_dataset()
ds_rlx[rlx_name].attrs['units'] = 's-1'
ds_rlx[rlx_name].attrs['cell_method'] = 'time: point'

cwd = os.getcwd()
ds_rlx.attrs["history"] = f"Created {cwd}/{btx}"


if f_save:
  encoding = {rlx_name: {'_FillValue': None}}
  flout = f'relax_rate_{int(rate_max_hrs):03d}hrs.nc'
  pthsis = gridfls['MOM6_NEP'][run_name]['pthsis']
  dflout = os.path.join(pthsis, flout)

  print(f'Saving SIS2 relaxation time scale --> {dflout}')
  ds_rlx.to_netcdf(
      dflout,
      format='NETCDF3_64BIT',
      engine='netcdf4',
      encoding=encoding
  )


if check_rlx:
  plt.ion()

  clrmp = mclrmps.colormap_temp2()
  clrmp = mclrmps.colormap_conc() 
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  cff = 1.e5
  

  # Stereographic Map projection:
  from mpl_toolkits.basemap import Basemap, cm
  m = Basemap(width=5000*1.e3,height=5000*1.e3, resolution='l',\
              projection='stere', lat_ts=50, lat_0=62, lon_0=-165)

  xR, yR = m(hlon, hlat)

  AP = RLXIS.copy()*cff
  AP = np.where(HH>=0., np.nan, AP)


  rmin = 0.
  rmax = 25.
  apmax = np.nanmax(AP)
  if apmax > 1.:
    rmax = np.floor(apmax)
  else:
    rmax = apmax 
  
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax0 = plt.axes([0.1, 0.1, 0.8, 0.8])
  m.drawcoastlines()
  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

  img = ax0.pcolormesh(xR, yR, AP, cmap=clrmp, vmin=rmin, vmax=rmax)
#  img = ax0.pcolormesh(RLXHR, cmap=clrmp)
  if rate_max_hrs <= 2:
    tscntrs = [1,2,5,10,40,80,120,240,360,480,600,720,1440]
  elif rate_max_hrs <=4:
    tscntrs = [4,6,10,40,60,80,120,240,360,480,600,720,1440]
  elif rate_max_hrs <=24:
    tscntrs = [24,26,30,60,80,120,240,360,480,600,720,960,1440]
  elif rate_max_hrs <=120:
    tscntrs = [120,150,240,360,480,600,720,960,1440]
  else:
    tscntrs = [120,150,240,360,480,600,720,960,1440]


  tslabels = tscntrs
  CS = ax0.contour(xR,yR,RLXHR,tscntrs, linestyles='solid', linewidths=1, colors=[(0., 0., 0.)])
  ax0.clabel(CS, tslabels,inline=1, fontsize=10)

  ax0.set_title(f'Relaxation rate (s-1), contours: hrs, strongest rlx {rate_max_hrs:.1f} hrs')

  ax2 = fig1.add_axes([ax0.get_position().x1+0.025, ax0.get_position().y0,
                     0.02, ax0.get_position().height])
  # extend: min, max, both
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
  ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)
  ax2.set_ylabel(f'Relaxation, {1./cff:.1e}, s-1')

  bottom_text(btx, pos=[0.2, 0.01])

# Plot mapping/transformations of the reference domain D 
# in order to prepare relaxation rate fields for the NEP domain
def axes_refdom(ax0,X,Y):
  ax0.axis('scaled')
  ax0.set_xlim([-np.max(X), np.max(X)])
  ax0.set_ylim([Y[0],Y[-1]])
  ax0.grid('on')

  return ax0

if check_ref_domain:
  """
   Plot remapping stages of the reference domain 
  """
  plt.ion()

  #clrmp = mclrmps.colormap_temp2()
  #rmin = 0.
  #rmax = 1.
  clrmp = mclrmps.colormap_conc()
  clrmp.set_bad(color=[1, 1, 1])

  fig2 = plt.figure(2,figsize=(9,8))
  plt.clf()
  ax21 = plt.axes([0.05, 0.55, 0.4, 0.4])
  img = ax21.pcolormesh(X,Y,RLX, cmap=clrmp)
  ax21 = axes_refdom(ax21,X,Y)
  ax21.set_title('Reference domain')

  
  ax22 = plt.axes([0.5, 0.55, 0.4, 0.4])
  ax22.pcolormesh(X,Y,RMAP1, cmap=clrmp)
  ax22 = axes_refdom(ax22,X,Y)
  ax22.set_title(f'Mapping 1: f(z)=z^(1/{rdnm})')

  ax23 = plt.axes([0.05, 0.05, 0.4, 0.4])
  ax23.pcolormesh(XR,YR,RMAP2, cmap=clrmp)
  ax23 = axes_refdom(ax23,X,YR)
  ax23.set_title('Mapping2: f(z)=z*exp(tht)')

  ax24 = plt.axes([0.55, 0.05, 0.02, 0.4])
  clb = plt.colorbar(img, cax=ax24, orientation='vertical', extend='both')
  #ax22.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  #ax22.set_yticklabels(ax22.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  #clb.ax.set_yticklabels(["{:.1f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)
  ax24.set_title('Relaxation rate, s-1') 
 
  btx = 'relax_timescale.py'
  bottom_text(btx, pos=[0.2, 0.01])


