"""
  Plot sea ice conc/thickness  
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
from copy import copy
import matplotlib.colors as colors
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

from mod_utils_fig import bottom_text
import mod_plot_xsections as mxsct
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
#import mod_valid_utils as mvutil
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_anls_seas as manseas
import mod_utils_ob as mutob
importlib.reload(mutob)

# experiment: year start, month start, ...
# change dayrun to plot desired date output - # of days since start date
# in daily-mean output fields: date is in the middle of the averaging period
varnm  = 'ithck'  # iconc or ithck 
f_cntrobs = True
f_obsthck = True

if not varnm == 'iconc':
  f_cntrobs = False
if not varnm == 'ithck':
  f_obsthck = False


# Start of the run - needed only for seasonal forecasts:
YRS    = 1993 # year start of the forecast
MOS    = 4
DDS    = 1    
nens   = 2    # ens # for ensemble runs
dnmbR  = mtime.datenum([1994,3,1])  # day to plot

expt    = "seasonal_fcst"
runname = f'NEPphys_frcst_climOB_{YRS}-{MOS:02d}-e{nens:02d}'
#expt    = 'NEP_BGCphys_GOFS'
#runname = 'NEP_physics_GOFS-IC'
#expt    = 'NEP_seasfcst_LZRESCALE'
#runname = 'NEPphys_LZRESCALE_climOB_1993_04-e02'
dnmbS   = mtime.datenum([YRS,MOS,DDS]) 
dvR     = mtime.datevec(dnmbR)
dnmb0   = dnmbR
dv0     = mtime.datevec(dnmb0)
YR0, MM0, DD0 = dv0[:3]
jday0   = int(mtime.date2jday([YR0,MM0,DD0]))

# For thickness - show observed values if available:
if f_obsthck and (YR0 == 1993 or YR0 == 1994) and ( MM0>=3 and MM0<=5):
  f_obsthck = True
else:
  f_obsthck = False
  

print(f'Expt: {expt} Run: {runname} Plot date: {dvR[0]}/{dvR[1]}/{dvR[2]}')

if expt == 'NEP_BGCphys_GOFS':
  outfld = 'ice'
else:
  outfld = 'icem'

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

if expt == 'seasonal_fcst':
  pthfcst = pthseas['MOM6_NEP'][expt]['pthoutp'].format(runname=runname)
  pthfcst = os.path.join(pthfcst,f'{outfld}_{dvR[0]}{dvR[1]:02d}')
else:
  pthfcst  = pthseas['MOM6_NEP'][expt]['pthoutp'].format(YY=YR0, MM=MM0)
pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)
ndav       = pthseas['MOM6_NEP'][expt]['ndav']  # # of days output averaged

# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

# Find closest output:
YR0, jday0, dnmb0, flname_out = manseas.find_closest_output(pthfcst, dnmbR, fld=outfld)
dv0  = mtime.datevec(dnmb0)
YR0, MM0, DD0 = dv0[:3]

  
flice_name = pthseas['MOM6_NEP'][expt]['ficename'].format(YR=YR0, jday=jday0)
dfsis2 = os.path.join(pthfcst, flice_name)

# Averaging period:
dnmb_av1 = dnmb0 - np.floor(ndav/2)
#if dnmb_av1 < dnmbS: dnmb_av1=dnmbS
dnmb_av2 = dnmb_av1 + ndav-1

dset   = xarray.open_dataset(dfsis2)

if varnm == 'iconc':
  A2d = dset['siconc'].isel(time=0).data
elif varnm == 'ithck':
  A2d = dset['sithick'].isel(time=0).data

# Add observed ice conc:
# fields are downloaded from the Near-Real-Time NOAA/NSIDC 
# Climate Data Record of Passive Microwave Sea Ice Concentration 
# https://nsidc.org/data/g10016
if f_cntrobs:
  pthnsidc = pthseas["NRT_NSIDC"]['pthdaily'].format(YR=YR0)
  flnsidc  = f'seaice_conc_daily_nh_{YR0}{MM0:02d}{DD0:02d}_f11_v04r00.nc'
  drflnsidc = os.path.join(pthnsidc, flnsidc)
  print(f'Reading NTR NSIDC ice conc: {drflnsidc}')
  dset_nsidc = xarray.open_dataset(drflnsidc)

  ICnrt = dset_nsidc['nsidc_nt_seaice_conc'].data[0,:].squeeze()
  Xnrt  = dset_nsidc['xgrid'].data
  Ynrt  = dset_nsidc['ygrid'].data

  ICnrt = np.where(ICnrt>1., np.nan, ICnrt)

# Flip NSIDC grid:
  ICnrt = np.flipud(ICnrt)
  Ynrt  = np.flipud(Ynrt)

# Convert Polar Coordinates to Geostatic coordinates (lon/lat)
  import mod_misc1 as mmisc
  XX, YY = np.meshgrid(Xnrt, Ynrt, indexing='xy')
  LONnrt, LATnrt = mmisc.convert_polarXY_lonlat(XX,YY)

# Make lon 0, 360 to match NEP grid
  LONnrt = np.where(LONnrt<0, LONnrt+360., LONnrt)

# Get ice edge contour in the Bering Sea
# get rid of ice in unneeded part of the domain
  ICnrt[:,150:] = np.nan
  ICnrt[:200,:] = np.nan

  CNTR = manseas.derive_ice_contour(ICnrt, nmin=5)

def color_scatter(J, I, HIo, clrmp, rmin, rmax):
  """
  Select colors for scatter plot
  """
  Nclrs = clrmp.N
  cval  = np.linspace(rmin, rmax, Nclrs)
  vclrs = np.zeros((len(J),3))
  for ii in range(len(J)):
    j0 = J[ii]
    i0 = I[ii]
    hh = HIo[j0,i0]
    if hh > cval[-1]:
      vclrs[ii,:] = clrmp(Nclrs)[:3]
    else:
      iclr = max(np.where(cval <= hh)[0])
      vclrs[ii,:] = clrmp(iclr)[:3]

  return vclrs

# -------------------
#
# Plot ice fields
#
# -------------------
if varnm == 'iconc':
  clrmp = mclrmps.colormap_conc()
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 0.
  rmax = 1.
elif varnm == 'ithck':
  clrmp = mclrmps.colormap_ice_thkn()
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 0.
  rmax = 5.

dv_av1 = mtime.datevec(dnmb_av1)
yrs, mms, dds = dv_av1[:3]
dv_av2 = mtime.datevec(dnmb_av2)
yre, mme, dde = dv_av2[:3]

sttl = f"{runname} {varnm} avrg: {yrs}/{mms}/{dds}-{yre}/{mme}/{dde}"
# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=3300*1.e3,height=3300*1.e3, resolution='l',\
            projection='stere', lat_ts=55, lat_0=62, lon_0=-175)

xR, yR = m(hlon, hlat)

# Observied ice thickness
if f_obsthck:
  if MM0 == 3:
# ERS1 mean March ice thickness 1993-2001
    pthdata = '/work/Dmitry.Dukhovskoy/data/NSIDC_icethckn'
    flithkn = 'ers1_seaice_thickness_mean_march_1993to2001.nc'
    dflithkn = os.path.join(pthdata, flithkn)
    varthck = 'thickness'
    sinfo2 = 'Color bullets:  ERS1 derived mean ice thickness March 1993-2001\n'     
  else:
# Submarine spring draft obs:
    pthdata = '/work/Dmitry.Dukhovskoy/data/NSIDC_icethckn'
    flithkn = f'sub_seaice_thickness_mean_spring_{YR0}.nc'
    dflithkn = os.path.join(pthdata, flithkn)
    varthck = 'thick'
    sinfo2 = f'Color bullets:  submarine mean ice thickness spring {YR0}\n'     
  dset = xarray.open_dataset(dflithkn)
  LAT  = dset['lat'].data
  LON  = dset['lon'].data
  HIo  = dset[varthck].data
  LON = np.where(LON<-900, np.nan, LON)
  LAT = np.where(LAT<-900, np.nan, LAT)
  LON  = np.where(LON<0., LON+360., LON)
  xRo, yRo = m(LON, LAT)
  VMsk = ( (xRo < 1.e20) & (yRo < 1.e20) )
  J, I = np.where( (~np.isnan(HIo))  &  (HIo>0.) )
  vclrs = color_scatter(J, I, HIo, clrmp, rmin, rmax)

sinfo = f'SIS2: {dfsis2}\n'
if f_obsthck:
  sinfo = sinfo + sinfo2 

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

img = m.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
#m.contour(xR, yR, HH, [-1000], colors=[(0,0,0)], linestyles='solid')
#ax1.axis('scaled')
ax1.set_title(sttl)

# Show observed ice edge:
import mod_rtofs as mrtofs
if f_cntrobs:
  ncc = len(CNTR)
  for icc in range(ncc):
    Ic = CNTR[icc][:,0]
    Jc = CNTR[icc][:,1]
# Get geodetic coordinates:
# use exact (float) indices to interpolate exact geodetic coordinate
    nic = len(Ic)
    Xc = np.zeros((nic))
    Yc = np.zeros((nic))
    for ipp in range(nic):
      ii0 = Ic[ipp]
      jj0 = Jc[ipp]
      xc0, yc0 = mrtofs.interp_indx2lonlat(ii0, jj0, LONnrt, LATnrt)
      Xc[ipp] = xc0
      Yc[ipp] = yc0

    Xcm, Ycm = m(Xc,Yc)
    ax1.plot(Xcm, Ycm, linewidth=2, color=[1,0,0])

  ss2   = f'NSIDC iconc: {drflnsidc}' 
  sinfo = sinfo + ss2

# Plot observed ice thickness
if f_obsthck:
  m.scatter(xRo[J,I], yRo[J,I], c=vclrs, s=22)


if varnm == 'ithck':
  ax1.contour(xR, yR, A2d,[6, 10, 14], linestyles='solid', colors=[(0.95, 0.95, 0.95)])

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
# extend: min, max, both
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.1f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

ax3 = fig1.add_axes([0.02, 0.025, 0.8, 0.05])
ax3.text(0, 0, sinfo, fontsize=8)
ax3.axis('off')


btx = 'plot_seaice.py'
bottom_text(btx, pos=[0.2, 0.01])


