"""
  Calculate modif. hausdorff distance statistics
  following Dukhovskoy et al., 2015
  for irlx expmt. vs PIOMAS

  ice edge contour

  NOAA NWS EMC Dmitry Dukhovskoy
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib
import xarray
import matplotlib.colors as colors
from yaml import safe_load
import argparse

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
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hausdorff')

import mod_time as mtime
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_sis2_relax as msisrlx
import mod_rtofs as mrtofs
import mod_hausdorff_distance as mmhd

parser = argparse.ArgumentParser()
parser.add_argument("--yr", help="year to plot, default 2001 for NEP and 1995 for ARC", type=int)
parser.add_argument("--regn", help="NEP or ARC", type=str, required=True)
parser.add_argument("--intrp", help=" =1: interp PIOMAS mnth to daily for better accur., default=1", \
                    type=int)
parser.add_argument(
    "--enmb",
    help="List of experiment numbers (e.g., 1 3 5 32 33)",
    type=int,
    nargs="+",             # <-- allows one or more integers and will generate a list
    required=True
)
args = parser.parse_args()

# Test runs were performed for only 1 year
regn = args.regn if args.regn else None
YRS = args.yr if args.yr else None
interp = args.intrp if args.intrp else 1
ENMBS = args.enmb if args.enmb else None

# Number of test runs:
Nexpts = len(ENMBS)

interp_mnthly = interp > 0  # for more accurate comparison, do time interpolation of PIOMAS 
                      # to get mnthly mean values, similar to how it is done in SIS2
                      # when deriving iconc ithkn for day=d0 from PIOMAS target fields

varnm = 'iconc'
if YRS is None:
  if regn == 'NEP':
    YRS = 2001
  else:
    YRS = 1995

mstart = 1  # current test runs all started on Jan 1, 2001
MMS = 1
MME = 12

#if expt_nmb is None:
#  EXPTS=[1,2,3,4,5]
#else:
#  EXPTS=[expt_nmb]

# relax hours:
#RLXH = [0,1,24,120,360]


fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

if regn == 'NEP':
  pthtopo    = pthseas['MOM6_NEP']['seasonal_daily']['pthgrid']
  fgrid      = pthseas['MOM6_NEP']['seasonal_daily']['fgrid']
  ftopo_mom  = pthseas["MOM6_NEP"]['seasonal_daily']["ftopo"]
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

else:
  pthdata = '/work/Dmitry.Dukhovskoy/ARC12/irlx'
  pthtopo = '/work/Dmitry.Dukhovskoy/ARC12/topo_grid'

  dflarc  = os.path.join(pthtopo,'ocean_hgrid.nc')
  hlon, hlat = mmom6.read_mom6grid(dflarc, grdpnt='hgrid')

  dtopo = os.path.join(pthtopo,'ocean_topog.nc')
  ds_topo = xarray.open_dataset(dtopo)
  HH = -ds_topo['depth'].data

jdm, idm = HH.shape

lon_offsetW = lon_offsetE = None
def region_lim(A2d, regn):
  if regn == 'NEP':
    lat_min = 56.5
    A2d[hlat <= lat_min] = np.nan
    A2d[:572,:] = np.nan
    A2d[:598,:180] = np.nan
    A2d[748:,:140] = np.nan 
    # Ignore values near the Nw and NE bndries:
    lon_offsetW = 164.
    lon_offsetE = 222.
  elif regn == 'ARC':
    lat_min = 56.
    #lat_min = 60.
    A2d[hlat <= lat_min] = np.nan
    A2d[:55,:] = np.nan
    A2d[172:286,396:] = np.nan
    A2d[544:,:95] = np.nan
    A2d[477:,:74] = np.nan
    #A2d[432:,351:] = np.nan 

  return A2d 

def indx2lonlat(CI):
  print('Mapping NEP ice cntr index ---> Lon/lat ...')
  XC, YC = [], []
  for ik in range(len(CI)):
    if ik%50 == 0:
      prcnt = float(ik)/float(len(CI))*100.
      print(f'  {prcnt:5.2f}% done ...')

    xrt, yrt = mrtofs.interp_indx2lonlat(CI[ik,0], CI[ik,1], hlon, hlat)
    XC.append(xrt)
    YC.append(yrt)

  XC = np.array(XC)
  YC = np.array(YC)
      
  return XC, YC


MHD = np.zeros((12,Nexpts))   # months x N expts
npnts_min = 1  # min # of points in the contour to keep
TM = []
Nrec = 0
YRA = YRS
Mold = -1
imo = -1
for MMA in range(MMS,MME+1):
  print(f"Processing {YRA}/{MMA} IceConc ...")

  dnmb0 = mtime.datenum([YRA, MMA, 15])
  imo += 1

  iexp = 0
  for kexpt in range(Nexpts):
    expt_nmb = ENMBS[kexpt]
    if regn == 'NEP':
      pthtest = f'/archive/Dmitry.Dukhovskoy/fre/NEP/test_ice_relax/NEPphys_expt{expt_nmb:02d}/{YRS}-{mstart:02d}'
      pthrlx  = '/work/Dmitry.Dukhovskoy/NEP_input/SIS2_relax'
    elif regn == 'ARC':
      pthtest = f'/archive/Dmitry.Dukhovskoy/fre/ARC12/test_ice_relax/ARCphys_expt{expt_nmb:02d}/{YRS}-{mstart:02d}'
      pthrlx  = '/work/Dmitry.Dukhovskoy/ARC12/irlx'
      npnts_min=20

    dcice = os.path.join(pthtest,f'ice_month.nc')
    print(f'MM={MMA}, Reading {dcice}')

    # Note assumed start month = 1, if not - update imo 
    imo = MMA-1

    with xarray.open_dataset(dcice) as ds:
      CInep = ds['siconc'].isel(time=imo).data.squeeze()

    CInep[HH>=-0.1]=np.nan
    if lon_offsetW is not None and lon_offsetE is not None:
      CInep = np.where(hlon < lon_offsetW, np.nan, CInep)
      CInep = np.where(hlon > lon_offsetE, np.nan, CInep)
    CInep = region_lim(CInep, regn)


    if regn == 'NEP':
      TCNT = msisrlx.derive_iconc_contour(CInep, ic0=0.15, npmin=npnts_min)
    elif regn == 'ARC':
      TCNT = msisrlx.derive_iconc_contour_ARC(CInep, ic0=0.15, npmin=npnts_min)

    XCnep, YCnep = indx2lonlat(TCNT)

    # Read PIOMAS 
    if MMA != Mold:
      YR1 = YRA
      YR2 = YR1+1

      if regn == 'NEP':
        pthrlx  = '/work/Dmitry.Dukhovskoy/NEP_input/SIS2_relax'
        flout = f'PIOMASv21_ithkn_iconc_{YR1}_{YR2}_monthly.nc'
      elif regn == 'ARC':
        pthrlx = '/work/Dmitry.Dukhovskoy/ARC12/irlx'
        flout = f'PIOMASv21_ARC12_ithkn_iconc_{YR1}_{YR2}_monthly.nc'

      dfpiomas = os.path.join(pthrlx, flout)

      if interp_mnthly:
        CIpms = msisrlx.mnthly_PIOMAS_linear_daily(dfpiomas,dnmb0,'iconc')
        #HIpms = msisrlx.mnthly_PIOMAS_linear_daily(dfpiomas,dnmb0,'ithkn')
      else:
        CIpms,_ = msisrlx.read_relax_piomas(dnmb0, pthrlx, 'iconc')
        #HIpms,_ = msisrlx.read_relax_piomas(dnmb0, pthrlx, 'ithkn')

      CIpms[HH>=-0.1]=np.nan
      if lon_offsetW is not None and lon_offsetE is not None:
        CIpms = np.where(hlon < lon_offsetW, np.nan, CIpms)
        CIpms = np.where(hlon > lon_offsetE, np.nan, CIpms)
      CIpms = region_lim(CIpms, regn)

      if regn == 'NEP':
        TCNTp = msisrlx.derive_iconc_contour(CIpms, ic0=0.15, npmin=npnts_min)
      elif regn == 'ARC':
        TCNTp = msisrlx.derive_iconc_contour_ARC(CIpms, ic0=0.15, npmin=npnts_min)

      XCpms, YCpms = indx2lonlat(TCNTp)

      Mold = MMA

    P = np.column_stack((XCnep,YCnep))
    Q = np.column_stack((XCpms,YCpms))
    mhdGS = mmhd.modifHD(P, Q, geo2cart=True)

    MHD[imo,iexp] = mhdGS
    print(f'MHD = {mhdGS:.2f}')
    iexp += 1 


clrmp = mclrmps.colormap_conc()
clrmp.set_bad(color=[0.2, 0.2, 0.2])
rmin = 0.
rmax = 1.
cntr_clr = [0.9, 0.9, 0.9]
hcntrs = [0.15]

f_check = False
if f_check:
  plt.ion()
  fig1 = plt.figure(1,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  ax1.contour(hlon, hlat, HH, [0], colors=[(0,0,0)])
  ax1.axis('scaled')
  ax1.set_ylim([50, 80])
  ax1.set_xlim([157, 235])
  ax1.plot(XCnep,YCnep,'.')
  ax1.plot(XCpms,YCpms,'.')


  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  img = ax1.pcolormesh(CInep, cmap=clrmp, vmin=rmin, vmax=rmax)
  CS = ax1.contour(CInep, hcntrs, linestyles='solid', colors=[cntr_clr], linewidths=1)
  ax1.contour(CIpms, hcntrs, linestyles='solid', colors=[(0,0.8,0.6)], linewidths=1)
  ax1.plot(TCNT[:,0],TCNT[:,1],'.')
  ax1.plot(TCNTp[:,0],TCNTp[:,1],'.')

  # For NEP:
  if regn == 'NEP':
    xl1 = 24
    xl2 = 342
    yl1 = 565
    yl2 = 816
    ax1.set_xlim([xl1,xl2])
    ax1.set_ylim([yl1,yl2])


  # extend: min, max, both
  ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                     0.02, ax1.get_position().height])
  if fgnmb>1:
    clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='max')
  else:
    clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')

  ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  btx = 'hausdorff_iconc_irlxtest.py'
  bottom_text(btx, pos=[0.2, 0.05])


# Plot MHD:
#CLR = [[0.,0.4,0.9],
#         [0.9,0.5,0],
#         [0.,0.9,0.7],
#         [1.,0.9,0],
#         [0.8,0.,0.5],
#         [0.7, 1, 0.2]]

ECOLR = msisrlx.irlx_tests_colors()

TMyr = [x for x in range(1,13)]

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.5, 0.8, 0.4])
hndls = []
lnw0 = 2
for kk in range(Nexpts):
  expt_nmb = ENMBS[kk]
  mhd_exp = MHD[:,kk]
  expt_name = msisrlx.irlx_tests_name(expt_nmb, regn)

  if expt_nmb == 1:
    # Control run
    marker = 'o'
    mksz = 4
  else:
    marker = None
    mksz = None

  clr = ECOLR[kk]
  ln1, = ax1.plot(TMyr, mhd_exp, 
         linestyle='-', linewidth=lnw0, marker=marker, markersize=mksz, color=clr, label=expt_name)
  hndls.append(ln1)


ax1.set_xticks(TMyr)
ax1.set_xlim([TMyr[0]-0.1,TMyr[-1]+0.1])
ax1.grid('on')
ax1.set_xlabel('Months')
ax1.set_ylabel('MHD score, km')
sttl = f'MHD ice edge btw irxl expt and PIOMAS, {regn}'
ax1.set_title(sttl)
 
ax2 = plt.axes([0.7, 0.25, 0.25, 0.18])
ax2.legend(handles=hndls, loc='upper right')
ax2.axis('off')




