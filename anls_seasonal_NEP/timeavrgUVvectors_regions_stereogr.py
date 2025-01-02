"""
  Plot 2D mean U/V vectors 
  stereographic projection
  in different regions: # CalCur - Calif Current region, Alaska - Alaska region, BerSea - Bering
  Following Stoke et al., 2015

  in the Calif. Current region, see analysis:
  Auad et al., 2011
  The California Current System in relation to the Northeast Pacific Ocean circulation
  https://www.sciencedirect.com/science/article/pii/S0079661111001157

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
importlib.reload(manseas)


# Initial date
# Look at ens run #1 - the only ens. that has 5-day av. output fields
expt     = 'seasonal_daily'  # seasonal forecasts with dailyOB from SPEAR
varnm    = 'u'  # temp (potential) / salin
#dnmbS    = mtime.datenum([2015,1,1])
# Averaging time period:
MMS   = 1    # f/cast init. month in each year, can be changed to months: 1, 4, 7, 10
YAVRG = [x for x in range(2011,2021)]
MAVRG = [1,2,3]  # months to average: Winter  JFM, Summer: JAS
#MAVRG = [7,8,9]  # months to average: Winter  JFM, Summer: JAS
regn_name = 'CalUnderCur' # CalCur - Calif Current region, Alaska - Alaska region, BerSea - Bering
                     # Following Stoke et al., 2015
                     # CalUnderCur - region ~35-45 N, 0-200 m 

# Average over a depth range:
#LRS = [1, 13]  # lr 11 = -25 m
#LRS = [21, 30]  # 50 - 150 m
#LRS = [30, 35]  # 93 - 153 m
LRS = [34, 36]  # 137 - 171 m  - for Under Current

nensR    = 1
expt_nmb = 2   # 2 - seas f/casts with dailyOB

expt_name = f'NEPphys_frcst_climOB{expt_nmb:02d}'
run_info = f'{expt_name} init MM={MMS} e{nensR:02d}, avrg {varnm}: {min(YAVRG)}-{max(YAVRG)} Mo: {min(MAVRG)}-{max(MAVRG)}'

print(f'Plotting {varnm} {expt_name} ')
print(f'{run_info}')

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

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

#mo_fcsts = manseas.yrmo_seasonal_fcst(YRS, MMS)

ocnfld = 'oceanm'
pthoutp0 = pthseas['MOM6_NEP'][expt]['pthoutp'].format(expt_nmb=expt_nmb)
YR=2011
MM=4
subdir=f'oceanm_{YR}{MM:02d}'
pthfcst0 = os.path.join(pthoutp0,f'{YR}-{MM:02d}-e01','history')
list_files = manseas.list_oceanice_files(pthfcst0, prefix=ocnfld, subdir=subdir)
pthfull = os.path.join(pthfcst0,subdir)
floceanm = list_files[0]
ZM = manseas.read_oceanm3D_field(pthfull, floceanm, 'zl', notime=False)
ZM = -abs(ZM)
lr1 = LRS[0]
lr2 = LRS[1]
zz1 = ZM[lr1-1]
zz2 = ZM[lr2-1]

def select_color(CLRS, UINT, uval):
  """
    Select color for given value
  """
  ncc  = CLRS.shape[0]
  nint = len(UINT)
  if nint != ncc+1:
    raise Exception(" N of colors does not match U intervals")

  if uval <= UINT[0]: return(CLRS[0])
  if uval >= UINT[-1]: return(CLRS[-1])

  dd   = UINT-uval
  ibin = min(np.argwhere(dd >= 0.))[0] - 1
  clr0 = list(CLRS[ibin,:])

  return clr0


Time = []
iyr  = 0
for YRS in (YAVRG):
  dnmbS    = mtime.datenum([YRS,MMS,1])
  pthfcst0 = os.path.join(pthoutp0,f'{YRS}-{MMS:02d}-e{nensR:02d}','history')
  UU, TM = manseas.monthly_depth_mean_from_Ndaily3D(pthfcst0, YRS, MMS, 'u', \
                                                 ocnfld, lr1, lr2, MAVRG=MAVRG)

  VV, _  = manseas.monthly_depth_mean_from_Ndaily3D(pthfcst0, YRS, MMS, 'v', \
                                                 ocnfld, lr1, lr2, MAVRG=MAVRG)

  if iyr == 0:
    U2d = UU.copy()
    V2d = VV.copy()
  else:
    U2d = U2d + UU
    V2d = V2d + VV

  Time = Time + TM

  iyr += 1

U2d = U2d/iyr
V2d = V2d/iyr

# Collocate U/V:
U2c = mmom6.collocateU2H(U2d, 'symmetr', f_land0 = False)
V2c = mmom6.collocateV2H(V2d, 'symmetr', f_land0 = False)

S2d = np.sqrt(U2c**2 + V2c**2)
U2c = np.where(HH>0., np.nan, U2c)
V2c = np.where(HH>0., np.nan, V2c)
S2d = np.where(HH>0., np.nan, S2d)

DV = mtime.datevec2D(Time)

II = pthseas['ANLS_NEP'][regn_name]['II']
JJ = pthseas['ANLS_NEP'][regn_name]['JJ']
xlim1 = min(II)
xlim2 = max(II)
ylim1 = min(JJ)
ylim2 = max(JJ)

# Plot boundaries of the region:
lon_s = hlon[ylim1, xlim1:xlim2+1]
lat_s = hlat[ylim1, xlim1:xlim2+1]

lon_n = hlon[ylim2, xlim1:xlim2+1]
lat_n = hlat[ylim2, xlim1:xlim2+1]

lon_w = hlon[ylim1:ylim2+1, xlim1]
lat_w = hlat[ylim1:ylim2+1, xlim1]

lon_e = hlon[ylim1:ylim2+1, xlim2]
lat_e = hlat[ylim1:ylim2+1, xlim2]

Xreg, Yreg = mmisc.connect_segments([lon_w, lon_n, lon_e, lon_s], \
                                    [lat_w, lat_n, lat_e, lat_s])

# Land mask:
Lmsk = np.where(HH<0.,1.,0.)

# Colormaps:
cmp_Lmsk = mclrmps.colormap_landmask()
cmps  = mclrmps.colormap_discrete()
CLRS  = cmps.colors
nclrs = CLRS.shape[0]

rmin = 0.
rmax = 0.1
dU   = (rmax-rmin)/nclrs
UINT = np.arange(rmin, rmax+dU, dU)

btx = 'timeavrgUVvectors_regions_stereogr.py'
sttl = f"{run_info} \n {regn_name}  z={zz1:7.1f}/{zz2:7.1f} m"

# Stereographic projection:
from mpl_toolkits.basemap import Basemap, cm
match regn_name:
  case 'CalCur':
    width  = 2800*1.e3
    height = 3600*1.e3
    lat0   = 33.
    lon0   = -123.
   # Plotting vectors:
    dltI = 200  # how far offshore (grid points) to show the vectors
    nx = 20     # # of vectors in X dir
    ny = 20     # # of vectors in Y dir
    uscale = 5000. 
  case 'CalUnderCur':
    width  = 800*1.e3
    height = 1400*1.e3
    lat0   = 40.
    lon0   = -124.5
   # Plotting vectors:
    dltI = 50  # how far offshore (grid points) to show the vectors
    nx = 20     # # of vectors in X dir
    ny = 20     # # of vectors in Y dir
    uscale = 2000. 

m = Basemap(width=width, height=height, resolution='l',\
            projection='stere', lat_ts=35, lat_0=lat0, lon_0=lon0)

xR, yR = m(hlon, hlat)
hcntrs = [x for x in range(-8000,-1,500)]

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()

ax1 = manseas.plot_stereogr_axis(fig1, m, xR, yR, Lmsk, cmp_Lmsk, rmin, rmax, \
                       btx=btx, sttl=sttl, clrbar=False)
plt.sca(ax1)
ax1.contour(xR, yR, HH, hcntrs, colors=[(0.8, 0.8, 0.8)], linestyles='solid', linewidths=1)

# Plot vectors:
# ==================
print('Plotting vectors')

ny = 20
Yvec = np.linspace(ylim1+10, ylim2-10, ny).astype(int)

# Find coast point along Yvec:
IIv = np.zeros((len(Yvec), nx)).astype(int)
JJv = np.zeros((len(Yvec), nx)).astype(int)
XXv = np.zeros((len(Yvec), nx))
YYv = np.zeros((len(Yvec), nx))
for ii in range(len(Yvec)):
  jj0 = Yvec[ii]
  dmm = HH[jj0,:]
  ii0 = min(np.where(dmm >= 0.)[0])
  xx = np.linspace(ii0-dltI, ii0-1, nx).astype(int)
  JJv[ii,:] = jj0
  IIv[ii,:] = xx
  XXv[ii,:] = hlon[jj0, xx]
  YYv[ii,:] = hlat[jj0, xx]
  
xp, yp = m(XXv, YYv)

import mod_draw_vector as mvec

vcol = [0., 0., 0.]
a1, a2 = XXv.shape
for jj in range(a1):
  for ii in range(a2):
    xs = XXv[jj,ii]
    ys = YYv[jj,ii]
    # Find degree/km for lat and lon:
    dlat = 1./(mmisc.dist_sphcrd(ys,xs,ys+1.,xs)*1.e-3)
    dlon = 1./(mmisc.dist_sphcrd(ys,xs,ys,xs+1)*1.e-3)

    i0 = IIv[jj,ii]
    j0 = JJv[jj,ii]
    uu = U2c[j0,i0]*uscale*dlon
    vv = V2c[j0,i0]*uscale*dlat
    xe = xs + uu
    ye = ys + vv 
    xvect = np.array([xs,xe])
    yvect = np.array([ys,ye])

    beam1, beam2 = mvec.arrow_vertices([xs,ys],[xe,ye])
    xb1 = beam1[:,0]
    yb1 = beam1[:,1]
    xb2 = beam2[:,0]
    yb2 = beam2[:,1]

    xvp, yvp   = m(xvect, yvect)
    xb1p, yb1p = m(xb1, yb1)
    xb2p, yb2p = m(xb2, yb2)

    sm = np.sqrt(U2c[j0,i0]**2 + V2c[j0,i0]**2)
    if sm <= 1.e-5: continue
    if np.isnan(sm): continue
    vcol = select_color(CLRS, UINT, sm)

    m.plot(xvp,  yvp,  linewidth=1, color=vcol)
    m.plot(xb1p, yb1p, linewidth=1, color=vcol)
    m.plot(xb2p, yb2p, linewidth=1, color=vcol)

# Plot NEP domain:
#xdom, ydom = m(Xreg, Yreg)
#m.plot(xdom, ydom, 'w-')


# Colorbar for discrete colors:
sm = matplotlib.cm.ScalarMappable(cmap=cmps)
ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
clb = plt.colorbar(sm, cax=ax2, orientation='vertical', extend='max')

ax2.yaxis.set_ticks(list(np.linspace(0,1,nclrs+1)))
clb.ax.tick_params(direction='in', length=12)

uvals = np.linspace(rmin, rmax, nclrs+1)
tick_lbl = []
for uvl in uvals:
  tick_lbl.append(f'{uvl:0.2f}')
clb.ax.set_yticklabels(tick_lbl)
#ax2.set_yticklabels(ax2.get_yticks())
#ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
#clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)



