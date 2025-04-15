"""
  Plot monthly mean Bering Strait depth-integrated U
  Extracted in calc_BerStrFlux.py 
  Saved indiviual years in separate files

  usage: run BerStr_UV.py --expt 3 --MMI 7 --YRS 1993 --YRE 2008 --MMS 1 --MME 3

  MMS can be > MME, e.g. for winter average: MMS=12, MME=2 to include Dec - Feb. 
  However, for MMI=1, Jan, Feb and Dec are not in cnosequtive order being 11 months apart
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import pickle
import xarray
from yaml import safe_load
import argparse
from matplotlib.patches import Polygon

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

parser = argparse.ArgumentParser()
parser.add_argument("--expt", help="experiment number: 1, ...", type=int)
parser.add_argument("--MMI", help="init month of f/cast, 1,4,7,10", type=int)
parser.add_argument("--YRS", help="start year of f/casts to derive cold poot stat, 1993, ...", type=int)
parser.add_argument("--YRE", help="end year of f/casts to derive cold pool stat, 1993, ...", type=int)
parser.add_argument("--MMS", help="start calendar month for averaging, 1,...,12", type=int)
parser.add_argument("--MME", help="end calendar month for averaging, >=MMS 1, ..., 12", type=int)
args = parser.parse_args()

# Default values: 
YRS    = 1993 # year start of the forecast
YRE    = 2008
MMI    = 4
MMS    = 1
MME    = MMS   # one month
nens   = 1    # ens # for ensemble runs - =1 for 3D ocean fields
expt_nmb = 3  # =2 - seas. f/casts no irelax, =3 - seas. f/casts with ice relax
Tref   = 0.   # Ref t for computing heat flux
Sref   = 34.8 # Ref S for FW flux

if args.expt:
  expt_nmb = args.expt
if args.YRS:
  YRS = args.YRS
if args.YRE:
  YRE = args.YRE
if args.MMI:
  MMI = args.MMI
if args.MMS:
  MMS = args.MMS
  MME = MMS
if args.MME:
  MME = args.MME

# Flags:
plot_obs = True # add obs.-based estimates from Woodgate 2018

expt_name = "seasonal_daily"
runname   = f"NEPphys_frcst_dailyOB-expt{expt_nmb:02d}"

if MMI==1 and YRS==1993:
  print(f"first start year for init month {MMI} is 1994, changing {YRS} to 1994")
  YRS = 1994

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

pthtopo    = pthseas['MOM6_NEP'][expt_name]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt_name]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt_name]["ftopo"]
hgrid      = xarray.open_dataset(os.path.join(pthtopo,fgrid))
hmask      = xarray.open_dataset(os.path.join(pthtopo, 'ocean_mask.nc'))
dstopo_nep = xarray.open_dataset(os.path.join(pthtopo, ftopo_mom))
dfgrid_mom = os.path.join(pthtopo, fgrid)

# Hgrid lon. lat:
hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')

HH = dstopo_nep['depth'].data
HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)

jdm, idm = HH.shape


# Define Bering Strait:
IIb = [201, 201]
JJb = [680, 667]   # <--- this will result in positive flux to the Arctic Oc. 
#JJb = [667, 680]  # <--- this will result in positive flux to the Ber. Sea

# Define segment lengths, norms, etc:
# positive: norm vector is to the left as follow the section
# negative: to the right
import mod_mom6_valid as mom6vld
DX, DY = mmom6.dx_dy(hlon,hlat)
IJ     = np.column_stack((IIb, JJb))
SGMT   = mmisc.define_segments(IJ, DX, DY, curve_ornt='positive', check_pole=False)
II     = SGMT.I_indx
JJ     = SGMT.J_indx

# Arrange time series in calendar months from f/cast months:
Tfcst = np.arange(MMI,MMI+12)
Tcal = np.mod(Tfcst,12)
Tcal = np.where(Tcal==0,12,Tcal)

if MMS <= MME:
  mav = np.arange(MMS,MME+1)  # months for seasonal averaging
else:
  mav = np.arange(1,MME+1)
  mav = np.append(mav,np.arange(MMS,13))
  
#iorder = Tcal
#ip = np.where(Tcal==12)[0][0]+1
#if ip < 12:
#  icol = np.arange(12)
#  iorder = icol+ip
#  iorder = np.where(iorder>=12, iorder-12, iorder)

Usum = None
pthanls = pthseas['MOM6_NEP'][expt_name]['pthanls'].format(expt_nmb=expt_nmb)
nyrs = YRE-YRS+1
icc = 0
for YRI in range(YRS,YRE+1):
  dflnm = os.path.join(pthanls,f'mnthly_BerSea_Fluxes_expt{expt_nmb:02d}_{YRI}{MMI:02d}.pkl')
  # Saved fields: VFlx,FWFlx,HFlx,UV2d,TT,SS,ZM,ZZ,Hbtm,LSgm,XX,YY
  print(f'Loading Fluxes {dflnm}')
  with open(dflnm, 'rb') as fid:
    _,_,_,U2d,_,_,ZM,ZZ,Hbtm,LSgm,XX,YY = pickle.load(fid)

  #U2d = np.expand_dims(U2d, axis=0)
  # Select needed months
  for imm in range(len(mav)):
    mo = mav[imm]
    jmo = np.where(Tcal == mo)[0][0]
    if Usum is None:
      Usum = U2d[jmo,:,:]
      Usum = np.expand_dims(Usum, axis=0)
    else:
      dmm =  U2d[jmo,:,:]
      dmm = np.expand_dims(dmm, axis=0)
      Usum = np.append(Usum, dmm, axis=0) 

  icc += 1


Umn = np.mean(Usum, axis=0)
Ustd = np.std(Usum, axis=0)
nrec = Usum.shape[0]

# Depth-integrated transport:
nlev, nsgm = Umn.shape
dZ = abs(np.diff(ZZ))
Trnsp = np.zeros((nsgm))
Tstd  = np.zeros((nsgm))
Tmonth = np.zeros((nrec,nsgm))
for kk in range(nsgm):
  aa = np.sum(Umn[:,kk]*dZ*LSgm[kk])
  Trnsp[kk] = aa*1e-6    # Sv

  # Find StDev:
  for irec in range(nrec):
    bb = np.sum(Usum[irec,:,kk]*dZ*LSgm[kk])
    Tmonth[irec,kk] = bb*1.e-6  # Sv

Tstd = np.nanstd(Tmonth, axis=0)
  
  

def find_angle_XEast(m,XXm,YYm,II,JJ):
  vi_x = XXm[-1]-XXm[0]
  vi_y = YYm[-1]-YYm[0]
  vi = np.sqrt(vi_x**2 + vi_y**2)
  vi_x = vi_x/vi
  vi_y = abs(vi_y/vi)  # point in northward dir.

  # find unit vector in X dir on MOM6:
  j0=JJ[1]; i0=II[1]
  pi1 = m(hlon[j0,i0+1],hlat[j0,i0+1])
  pi_x = pi1[0]-XXm[0]
  pi_y = -(pi1[1]-YYm[0])
  pi = np.sqrt(pi_x**2 + pi_y**2)
  pi_x = pi_x/pi
  pi_y = pi_y/pi

  # Unit vectors, |a|=|b|=1
  aa = np.array([vi_x,vi_y])
  bb = np.array([pi_x,pi_y])
  cos_alfa = np.dot(aa,bb)
  alfa = np.arccos(cos_alfa)

  return alfa

def vector_stere(m,uu,vv,xs,ys,cff):
  """
   Plot a vector in stereographic projection
   The mapping tool does not preserve vector magnitudes
   when plotting the same vector pointing at different directions
  """
  # First find lengths of the vector pointing X or Y
  U = np.sqrt(uu**2 + vv**2)
  xe = xs + U*cff
  ye = ys  
  xvect = np.array([xs,xe])
  yvect = np.array([ys,ye])
  beam1, beam2 = mvec.arrow_vertices([xs,ys],[xe,ye], beta=20, larr_max=0.06, larr_min=0.02)
  xb1 = beam1[:,0]
  yb1 = beam1[:,1]
  xb2 = beam2[:,0]
  yb2 = beam2[:,1]
  xvp, yvp   = m(xvect, yvect)
  xb1p, yb1p = m(xb1, yb1)
  xb2p, yb2p = m(xb2, yb2)
  Lvec0  = np.sqrt(np.diff(xvp)**2 + np.diff(yvp)**2)[0]
  Lbeam0 = np.sqrt(np.diff(xb1p)**2 + np.diff(yb1p)**2)[0] 
  xv0 = np.diff(xvp)[0]
  yv0 = np.diff(yvp)[0]

  # Get projection of the true vector:
  xe = xs + uu*cff
  ye = ys + vv*cff
  xvect = np.array([xs,xe])
  yvect = np.array([ys,ye])
  beam1, beam2 = mvec.arrow_vertices([xs,ys],[xe,ye], beta=20, larr_max=0.06, larr_min=0.02)
  xb1 = beam1[:,0]
  yb1 = beam1[:,1]
  xb2 = beam2[:,0]
  yb2 = beam2[:,1]
  xvp, yvp   = m(xvect, yvect)
  xb1p, yb1p = m(xb1, yb1)
  xb2p, yb2p = m(xb2, yb2)
  Lvec   = np.sqrt(np.diff(xvp)**2 + np.diff(yvp)**2)[0]
  Lbeam1 = np.sqrt(np.diff(xb1p)**2 + np.diff(yb1p)**2)[0]
  Lbeam2 = np.sqrt(np.diff(xb2p)**2 + np.diff(yb2p)**2)[0]
  
  crct_cf = Lvec0/Lvec
  xvp[1] = xvp[0]+crct_cf*(xvp[1]-xvp[0])
  yvp[1] = yvp[0]+crct_cf*(yvp[1]-yvp[0])

  cbm1 = Lbeam0/Lbeam1
  cbm2 = Lbeam0/Lbeam2
  dx1 = xb1p[1]-xb1p[0]
  dy1 = yb1p[1]-yb1p[0]
  dx2 = xb2p[1]-xb2p[0]
  dy2 = yb2p[1]-yb2p[0]
  xb1p[1] = xb2p[1] = xvp[1]
  yb1p[1] = yb2p[1] = yvp[1]
  xb1p[0] = xb1p[1] - cbm1*(dx1)
  yb1p[0] = yb1p[1] - cbm1*(dy1)
  xb2p[0] = xb2p[1] - cbm2*(dx2)
  yb2p[0] = yb2p[1] - cbm2*(dy2)

  return xvp, yvp, xb1p, yb1p, xb2p, yb2p

# Plot vertical section of U northward:
from matplotlib.patches import Polygon
plt.ion()

btx = 'BerStr_UV.py'

clrmp_dlt = mclrmps.colormap_ssh(cpos='YlOrRd', cneg='PuBuGn_r')
clrmp_dlt.set_bad(color=[0.6,0.6,0.6])
rmin = -0.5
rmax = 0.5

scntrs = [x/100 for x in range(5,50,5)]

# Patch bottom:
verts = [(np.max(XX),-6000),*zip(XX,Hbtm),(np.min(XX),-6000)]
poly = Polygon(verts, facecolor='0.3', edgecolor='0.3', zorder=2)

sttl = f'Bering Str Northrward U & StDev\n {runname} init M={MMI:02d} avrg: {YRS}-{YRE}, {MMS}/{MME}'
clr_obs=[0,0,0]

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.07, 0.5, 0.5, 0.4])
img = ax1.pcolormesh(XX,ZM, Umn, cmap=clrmp_dlt, vmin=rmin, vmax=rmax)
ax1.set_ylim([-60,0])
ax1.set_xlim([190.2, 192.1])
CS = ax1.contour(XX,ZM,Ustd,scntrs, linestyles='solid', colors=[(0.,0.,0)])
ax1.clabel(CS, inline=True, fontsize=10)
ax1.add_patch(poly)
ax1.set_title(sttl)

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
# extend: min, max, both
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.1f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

# Plot mean depth-intgr transports:
# Stereographic Map projection:
LMsk = np.where(HH>=0, 0, 1)
cmp_Lmsk = mclrmps.colormap_landmask(clr0=[0.7,0.7,0.7])
from mpl_toolkits.basemap import Basemap, cm
m = Basemap(width=100*1e3, height=100*1e3, resolution='l',\
            projection='stere', lat_ts=67, lat_0=66.2, lon_0=-168.)

ax3 = plt.axes([0.12,0.07,0.35,0.35])
ax3.cla()

xR, yR = m(hlon, hlat)
XXm, YYm = m(XX,YY)
ax3.pcolormesh(xR,yR,LMsk, cmap=cmp_Lmsk, vmin=0, vmax=0.1)
CS1 = ax3.contour(xR,yR,HH,[-100,-75,-50,-25], linestyles='solid', colors=[(0.8, 0.8, 0.8)])
ax3.clabel(CS1, inline=True, fontsize=8)
ax3.plot(XXm,YYm,'.-')
ax3.axis('scaled')
xlim=0.08e6
ax3.set_xlim([-xlim,xlim])
ax3.set_ylim([-xlim,xlim])
#m.drawcoastlines()
#m.drawparallels(np.arange(-90.,120.,10.))
#m.drawmeridians(np.arange(-180.,180.,10.))

# Assuming northward transport, i.e. only 1 component
# Find i,j unit vectors matching the along- and cross-strait directions
alfa = find_angle_XEast(m,XXm,YYm,II,JJ)

import mod_draw_vector as mvec
cff = 3.5
vcol = [0., 0., 0.]
for isgm in range(nsgm):
  Utr = Trnsp[isgm]  # along X-axis in MOM6
  Vtr = 0.
  if abs(np.sqrt(Utr**2+Vtr**2)) < 1.e-3:
    continue
  # Project on local N/E dir:
  uu = Utr*np.cos(alfa) - Vtr*np.sin(alfa)
  vv = Utr*np.sin(alfa) + Vtr*np.cos(alfa)
  xs = XX[isgm]
  ys = YY[isgm]
  xvp, yvp, xb1p, yb1p, xb2p, yb2p = vector_stere(m,uu,vv,xs,ys,cff) 
  m.plot(xvp,  yvp,  linewidth=2, color=vcol)
  m.plot(xb1p, yb1p, linewidth=2, color=vcol)
  m.plot(xb2p, yb2p, linewidth=2, color=vcol)

  # Show StDev:
  umin = Trnsp[isgm] - Tstd[isgm]
  vmin = 0. 
  uu = umin*np.cos(alfa) - vmin*np.sin(alfa)
  vv = umin*np.sin(alfa) + vmin*np.cos(alfa)
  xminm, yminm,_ ,_ ,_ ,_ = vector_stere(m,uu,vv,xs,ys,cff)
 
  umax = Trnsp[isgm] + Tstd[isgm]
  vmax = 0. 
  uu = umax*np.cos(alfa) - vmax*np.sin(alfa)
  vv = umax*np.sin(alfa) + vmax*np.cos(alfa)
  xmaxm, ymaxm,_ ,_ ,_ ,_ = vector_stere(m,uu,vv,xs,ys,cff)

  col_std = [0.9, 0.4,0]
  m.plot([xminm[1],xmaxm[1]],[yminm[1],ymaxm[1]],'-', linewidth=1, color=col_std)
  m.plot(xminm[1],yminm[1],'o',markersize=3, markerfacecolor=col_std, markeredgecolor='none')
  m.plot(xmaxm[1],ymaxm[1],'o',markersize=3, markerfacecolor=col_std, markeredgecolor='none')

# Draw a unit vector:
isgm = 3
uu = 0.1
vv = 0.
xs = 189.9
ys = 65.2
xvp, yvp, xb1p, yb1p, xb2p, yb2p = vector_stere(m,uu,vv,xs,ys,cff)
m.plot(xvp,  yvp,  linewidth=2, color=vcol)
m.plot(xb1p, yb1p, linewidth=2, color=vcol)
m.plot(xb2p, yb2p, linewidth=2, color=vcol)
xt1,yt1 = m(xs,ys-0.1)
ax3.text(xt1,yt1,f'{uu:.2f} Sv')

ax3.set_title('Depth-integrated transport & StDev, Sv')
ax3.set_xticks([])
ax3.set_yticks([])

ax3.set_xlim([-xlim,xlim])
ax3.set_ylim([-xlim,xlim])

bottom_text(btx, pos=[0.05, 0.02])






