"""
  Plot 2D fields from GLORYS monthly fields 
  for visual comparison with NEP nudging restart file
  how close the nudged simulations are to GLORYS

Data fields from Liz:
location for monthly GLORYS means for NEP region:
/archive/e1n/datasets/GLORYS/monthly_means/

location for padded monthly GLORYS means, concatenated by year used to generate NEP clim nudging files:
/archive/e1n/datasets/GLORYS/monthly_climatologies/

location for monthly GLORYS means, regridded to NEP for nudging as individual months:
/archive/e1n/mom6/NEP/sponge/monthly_sponge_files/

location for padded monthly GLORYS means, regridded to NEP for nudging and concatenated by year:
/archive/e1n/mom6/NEP/sponge/clims/

The last directory contains the files I used for nudging the solution to GLORYS. 

Daily GLORYS reanalysis for NEP domain prepared by Liz:
/archive/e1n/datasets/GLORYS/YYYY/nep_10

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
import pickle
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
sys.path.append(PPTHN + '/TEOS_10/gsw')
sys.path.append(PPTHN + '/TEOS_10/gsw/gibbs')
sys.path.append(PPTHN + '/TEOS_10/gsw/utilities')
import mod_swstate as msw
import conversions as gsw

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

dnmb = mtime.datenum([1994,4,1])
expt = 'GLORYS_NEP'  # GLORYS extracted for NEP domain
varnm = 'zos'  # thetao, so, zos 
lr0  = 1  # ocean layers from 1, ..., 50

dv0 = mtime.datevec(dnmb)
YR0, MM0, DD0 = dv0[:3]

print(f'Plotting GLORYS {varnm} {YR0}/{MM0}/{DD0}')

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

pthdata = gridfls[expt]["monthly"]['pthoutp']
flnm = gridfls[expt]["monthly"]['fdata'].format(year=YR0)
dfl_glorys = os.path.join(pthdata,flnm)

# GLORYS monthly fields interpolated onto MOM6 NEP grid, use MOM topo
pthtopo    = gridfls['MOM6_NEP']['seasonal_fcst']['pthgrid']
fgrid      = gridfls['MOM6_NEP']['seasonal_fcst']['fgrid']
ftopo_mom  = gridfls["MOM6_NEP"]["seasonal_fcst"]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)

HH = dstopo_nep['depth'].data
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)


# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

dset = xarray.open_dataset(dfl_glorys)

# Select time: data padded, i.e. start from 0 and + 1 month at the end:
#  time = 0, 15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 
#    349.5, 380.5 ;
itime = MM0
idepth = lr0-1

# Subtract mean:
if varnm == 'zos':
  A2d = dset[varnm].isel(time=itime).data.squeeze()
  zmin = -200.
  dmm  = A2d.copy()
  dmm  = np.where(HH>zmin, np.nan, dmm)
  A2d = A2d - np.nanmean(dmm)
else:
  A2d = dset[varnm].isel(time=itime, depth=idepth).data.squeeze()

#rmin, rmax = mutob.minmax_clrmap(A2d, cpnt=0)

if varnm == 'so':
  clrmp = mutil.colormap_salin(clr_ramp=[1,0.85,1])
  clrmp.set_bad(color=[0.2, 0.2, 0.2])
  rmin = 30.0
  rmax = 35.0
elif varnm == 'thetao':
  clrmp = mutil.colormap_temp(clr_ramp=[0.9,0.8,1])
  clrmp.set_bad(color=[1,1,1])
  rmin = -2.
  rmax = 23.
elif varnm == 'zos':
  clrmp = mutil.colormap_ssh(nclrs=200)
  rmin = -0.5
  rmax = 0.5


# Set up orthographic projection
from mpl_toolkits.basemap import Basemap, cm

# Add extra row/col for plotting
# Add extra row/col for plotting with pcolormesh
#lonw = LONG.copy()
#latw = LATG.copy()
#lonw = np.insert(lonw, -1, lonw[:,-1]+0.01, axis=1)
#lonw = np.insert(lonw, -1, lonw[-1,:]+0.01, axis=0)
#latw = np.insert(latw, -1, latw[-1,:]+0.01, axis=0)
#latw = np.insert(latw, -1, latw[:,-1]+0.01, axis=1)


lon0 = 220.
lat0 = 50.
res  = 'l'
m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution=res)
xR, yR = m(hlon, hlat)
PMsk = ( (xR > 1e20) | (yR > 1e20) )
AA = A2d.copy()
AA = np.where(HH >= 0, np.nan, AA)
#AA = np.insert(AA, 0, AA[:,-1], axis=1)
#AA = np.insert(AA, -1, AA[-1,:], axis=0)
AA[PMsk] = np.nan
xR[PMsk]   = 1.e30
yR[PMsk]   = 1.e30

ny, nx = A2d.shape
#AA = AA[0:ny, 0:nx]
#AA = np.where(HHG >= 0., np.nan, AA)
#xBND, yBND = m(IBND, JBND)

#import mod_colormaps as mclrmp

plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
m.drawcoastlines()
im1 = m.pcolormesh(xR, yR, AA, cmap=clrmp, vmin=rmin, vmax=rmax)

m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(-180.,180.,10.))

#m.plot(xBND, yBND, 'r.')

sttl = f"GLORYS12v1, Monthly mean {YR0}/{MM0} {varnm}"
ax1.set_title(sttl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
# extend: min, max, both
clb = plt.colorbar(im1, cax=ax2, orientation='vertical', extend='both')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

sinfo=dfl_glorys
ax3 = fig1.add_axes([0.1,0.05,0.8,0.02])
ax3.text(0,0, sinfo, fontsize=8)
ax3.axis('off')


btx = 'plot_GLORYS_TS_ortho.py'
bottom_text(btx, fsz=6, pos=[0.05, 0.03])














