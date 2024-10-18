"""
  Plot T/S diagrams from WOA23
  Bering Sea
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import netCDF4
from netCDF4 import Dataset as ncFile
import importlib
import yaml
from yaml import safe_load

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
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_anls_seas as manseas
import mod_plot_xsections as mxsct

YR      = 1993
MM      = 1  # for seasonal indicate month, for annual: MM = 13

grd=0.25
if grd==0.25:
  cgrd=4
woa='woa23'

seas, decade, yr1_dec, yr2_dec = manseas.season_decade_woa(YR,MM) 

woa_seas = {"13": "Jan-Mar",
            "14": "Apr-Jun",
            "15": "Jul-Spt",
            "16": "Oct-Dec",
            "0": "annual"}

seas_nm = woa_seas[f"{seas}"]
print(f'Plotting T/S diagr for WOA23 seas={seas} {seas_nm} decade={decade}')

urlBase = 'https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/'
urlT    = f"{urlBase}temperature/netcdf/{decade}/0.25/"
urlS    = f"{urlBase}salinity/netcdf/{decade}/0.25/"
tfnm    = f"woa23_{decade}_t{seas:02d}_{cgrd:02d}.nc"
sfnm    = f"woa23_{decade}_s{seas:02d}_{cgrd:02d}.nc"
# All period average:
#tfnm=f'{woa}_decav_t{seas:02d}_{cgrd:02d}.nc'
#sfnm=f'{woa}_decav_s{seas:02d}_{cgrd:02d}.nc'

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

def read_field(furl,varnm):
  print("Reading {1} from {0}".format(furl,varnm))
  nc=ncFile(furl)
# lookup a variable
  dmm0 = nc.variables[varnm][:].data.squeeze()
  dmm = np.copy(dmm0)
  return dmm

def lookup_ncvar(nc):
  ii=0
  for var in nc.variables.values():
    ii+=1
    print('--------\n')
    print('Var # {0}'.format(ii))
    print(var)
 
def read_TS(furl, var_read):
  A3d = read_field(furl,var_read)
  A3d = np.where(A3d > 1.e10, np.nan, A3d) 

  latW = read_field(furl,'lat')
  lonW = read_field(furl,'lon')

  # reshaffle to have -180/180 lon inside the domain
  A3d, lonW = mmisc.shuffle3D_lon180(A3d, lonW) 
  lonW = np.where(lonW<0., lonW+360., lonW)

  return A3d, lonW, latW

furl = os.path.join(urlT, tfnm)
T3d, lonW, latW = read_TS(furl, 't_an')
furl = os.path.join(urlS, sfnm)
S3d,  _, _  = read_TS(furl, 's_an')
ZZ  = read_field(furl,'depth')
ZZ  = -abs(ZZ)

LONW, LATW = np.meshgrid(lonW, latW)

# Bering Sea domain, use NEP indices
# Find corresponding WOA indices
# Get NEP indices of the polygon:
import mod_rtofs as mrtofs
II = pthseas['ANLS_NEP']['poly_BerSea']['II']
JJ = pthseas['ANLS_NEP']['poly_BerSea']['JJ']
# Hgrid lon. lat:
expt    = "seasonal_fcst"
pthtopo    = pthseas['MOM6_NEP'][expt]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt]['fgrid']
dfgrid_mom = os.path.join(pthtopo, fgrid)
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')
XMOM = hlon[JJ,II]
YMOM = hlat[JJ,II]
XMOM = np.where(XMOM < 0., XMOM+360., XMOM)
IW = []
JW = []
for ik in range(len(XMOM)):
  ii, jj = mrtofs.find_indx_lonlat(XMOM[ik], YMOM[ik], lonW, latW) 
  IW.append(ii)
  JW.append(jj)


kdm, jdm, idm = T3d.shape
# Construct approximate depths:
HHW = np.zeros((jdm,idm))-6000.  # WOA deepest depth is -5500
for kk in range(1,kdm):
  a1 = T3d[kk-1,:].squeeze()
  a2 = T3d[kk,:].squeeze()
  HHW = np.where( ((np.isnan(a2)) & (~np.isnan(a1))), ZZ[kk-1], HHW) 
a1 = T3d[0,:].squeeze()
HHW = np.where(np.isnan(a1), 10., HHW)

DX, DY = mmom6.dx_dy(LONW, LATW)
Acell  = DX*DY
X, Y   = np.meshgrid(np.arange(idm), np.arange(jdm))
MS, _, _ = mmisc.inpolygon_v2(X, Y, IW, JW)  # 
JBS, IBS = np.where( (MS == 1) & (HHW >= -250) & (HHW < 0) ) #exclude deeep regions
MSKBS  = np.zeros((jdm,idm))
MSKBS[JBS,IBS] = 1

# Convert in situ --> potential T at z=0:
import mod_regmom as mregmom
print('Converting T in situ --> T potential')
kdm  = len(ZZ)
Z3d  = np.tile(ZZ, idm*jdm).reshape((idm,jdm,kdm))
Z3d  = np.transpose(Z3d, (2, 1, 0))
Tp3d = mregmom.insitu2pot_3D(T3d, S3d, Z3d, LATW)

# Compute conservative T from potential T
#Tcons = gsw.CT_from_pt([20, 0, 35],[20, 0, 25]) 
# Derive potential enthalpy as Hpot = Tcons*Cp0 
# where Cp0 is the heat capacity of seawater to be 
# 53989.244 952 928 15 J kg-1 K-1
# potential enthalpy [J/kg] can be used as "heat content per unit mass"
# Compute absolute salinity from practical S:
sys.path.append(PPTHN + '/TEOS_10/gsw')
sys.path.append(PPTHN + '/TEOS_10/gsw/gibbs')
sys.path.append(PPTHN + '/TEOS_10/gsw/utilities')
import mod_swstate as msw
import conversions as gsw
PR   = np.zeros((kdm,jdm,idm))
for kk in range(kdm):
  pr_db, _ = msw.sw_press(Z3d[kk,:,:].squeeze(), LONW)
  PR[kk,:] = pr_db

# Estimate layer thickness:
dP = manseas.derive_dP_WOA(ZZ, Tp3d)

print('Computing absolute Salinity')
SA = gsw.SA_from_SP(S3d, PR, LONW, LATW)

# Discard regions outside Bering Sea:
# Fill bottom with nans
for ik in range(kdm):
  a2d = Tp3d[ik,:]
  a2d = np.where(MSKBS < 1, np.nan, a2d)
  a2d = np.where(HHW > ZZ[ik], np.nan, a2d)
  Tp3d[ik,:] = a2d

  s2d = SA[ik,:]
  s2d = np.where(MSKBS < 1, np.nan, s2d)
  s2d = np.where(HHW > ZZ[ik], np.nan, s2d)
  SA[ik,:] = s2d

  Zpnts = np.where(np.isnan(s2d), np.nan, ZZ[ik])/abs(HHW)  # normalize by depth
  Z3d[ik,:] = Zpnts

T1d = Tp3d[(~np.isnan(Tp3d))]
S1d = SA[(~np.isnan(SA))]
Z1d = Z3d[(~np.isnan(Z3d))]

# Density
SS, TT  = np.meshgrid(np.arange(15,38,0.1), np.arange(-1.8,12,0.1))
sgm0 = msw.sw_dens0(SS, TT)-1000.

#
# Colormap
import mod_colormaps as mclrmps
clrmp = mclrmps.colormap_temp(clr_ramp=[0.4,0.,0.6])
rmin  = -1.
rmax  = 0.

seas_nm = woa_seas[f"{seas}"]
sttl = f"WOA23 Bering Sea shelf,  decade:{yr1_dec}-{yr2_dec} {seas_nm}"
furl = os.path.join(urlT, tfnm)
ss1 = f"T: {furl}\n"
furl = os.path.join(urlS, sfnm)
ss2 = f"S: {furl}\n"
ss3 = "in situ T --> pot. T, Pract S --> abs. S"
sinfo = ss1 + ss2 + ss3


# Plot T/S
plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.2, 0.3, 0.6, 0.6])
img = ax1.scatter(S1d, T1d, c=Z1d,  vmin=rmin, vmax=rmax, s=8, cmap=clrmp)

ax1.axis('scaled')
ax1.set_xlim([28, 35.3])
ax1.set_ylim([-2, 11])

CS = ax1.contour(SS, TT, sgm0, [x for x in range(0,32,1)], linestyles='solid',
            colors=[(0.6,0.6,0.6)])
ax1.clabel(CS, inline=True, fontsize=11)

#ax1.set_xticks(Tticks)
#ax1.grid(True)
#ax1.set_xticklabels(tck_lbls)
ax1.set_ylabel(r'$\theta$$^oC$')
ax1.set_xlabel('Absolute S')
ax1.set_title(sttl)

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

ax3 = fig1.add_axes([0.02, 0.15, 0.8, 0.05])
ax3.text(0, 0, sinfo, fontsize=10)
ax3.axis('off')

btx = 'plot_TSdiagr_WOA23.py'
bottom_text(btx, pos=[0.1, 0.1])


