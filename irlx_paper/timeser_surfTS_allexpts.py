"""
  Show timeseries of monthly surface T/S
  from irlx experiments 
  to check for the drifts due to relaxation

  NOAA NWS EMC Dmitry Dukhovskoy
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import xarray
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
import mod_mom6 as mmom6
import mod_sis2_relax as msisrlx

parser = argparse.ArgumentParser()
parser.add_argument("--yr", help="year to plot, default 2001 for NEP and 1995 for ARC", type=int)
parser.add_argument("--regn", 
                   help="subregion to analyze: NEP10k Arctic, NEP10k Bering Sea,  ARC10k", 
                   choices=['NEPA','NEPB','ARC'],
                     type=str, required=True)
parser.add_argument("--varnm", help="temp or salin", type=str, required=True)
args = parser.parse_args()

# Test runs were performed for only 1 year
regn = args.regn if args.regn else None
YRS = args.yr if args.yr else None
varnm = args.varnm if args.varnm else None

if YRS is None:
  if regn == 'NEPA' or regn == 'NEPB':
    YRS = 2001
  elif regn == 'ARC':
    YRS = 1995

mstart = 1  # current test runs all started on Jan 1, 2001
MMS = 1
MME = 12
EXPTS=[1,2,3,4,5]

# relax hours:
RLXH = [0,1,24,120,360]


fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

if regn == 'NEPA' or regn == 'NEPB':
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

if regn == 'NEPA' or regn == 'NEPB':
  # Bering Sea - Chukchi Sea :
  hsh = -5000.
  LMsk = np.where((HH>=hsh) & (HH<0), 1, 0)
  # Mask out southern lats:
  LMsk = np.where(hlat<55.,0,LMsk)
  LMsk[:567,:] = 0
  LMsk[:,:39] = 0
  LMsk[:595,177:] = 0
  LMsk[:579,:129] = 0
  LMsk[:575,:143] = 0
  #LMsk[748:,:143] = 0

  # Remove near-boundary points:
  LMsk[810:,:] = 0
  LMsk[:,338:] = 0

  # 
  # Mask for Bering Sea
  # Bounded by the Bering Strait 
  BMsk = LMsk.copy()
  BMsk = np.where(hlat>66,0,BMsk)
  JB,IB = np.where(BMsk==1)
  # Mask for the Arctic Oc. part of the domain:
  # Ber. Str. + S. Chukchi Shelf
  AMsk = LMsk.copy()
  AMsk = np.where(BMsk==1, 0, AMsk)
  AMsk[:,:192] = 0


def region_lim(A2d, regn):
  if regn == 'NEPB':
   # Bering Sea / Chukchi sea regions in NEP10k
   # similar to statistcs analysis 
   # Use Bering Sea only
   A2d[BMsk==0] = np.nan     # Bering Sea only
  elif regn == 'NEPA':
   # Bering Sea / Chukchi sea regions in NEP10k
   # Use Chukchi Sea only
   A2d[AMsk==0] = np.nan    # average over the Arctic portion of the domain
  elif regn == 'ARC':
    lat_min = 60.
    A2d[hlat <= lat_min] = np.nan
    A2d[:55,:] = np.nan
    A2d[172:286,396:] = np.nan
    A2d[544:,:95] = np.nan
    A2d[477:,:74] = np.nan
    #A2d[432:,351:] = np.nan 

  return A2d 


DX, DY = mmom6.dx_dy(hlon, hlat)
Acell  = DX*DY*1.e-6  # km2

if regn == 'NEPA' or regn == 'NEPB':
  if varnm == 'temp':
    varnc = 'potT'
  elif varnm == 'salin':
    varnc = 'salt'
else:
  if varnm == 'temp':
    varnc = 'tos'
  elif varnm == 'salin':
    varnc = 'sos'

DAYS_OCN = [x for x in range(3,366,5)] # 5-day ocean output fields
nrecs = len(DAYS_OCN)

Acell_msk = None
nexpts = len(EXPTS)
TM = []
FLD = np.zeros((nrecs,nexpts))
irec = -1
for jday in DAYS_OCN:
  YR0 = YRS  # asuuming 1-yr run
  dnmb0 = mtime.jday2dnmb(YR0, jday)
  YR,MM,DD = mtime.datevec(dnmb0)[:3]
  assert YR==YR0, f'Check date conversion years do not match {YR0} and {YR}'
  print(f"Processing {YR}/{MM}/{DD} {varnm} ...")
  irec += 1

  TM.append(dnmb0)

  iexp = 0
  for expt_nmb in EXPTS:
    if regn == 'NEPA' or regn == 'NEPB':
      pthtest = f'/archive/Dmitry.Dukhovskoy/fre/NEP/test_ice_relax/NEPphys_expt{expt_nmb:02d}/{YRS}-{mstart:02d}'
      docn = os.path.join(pthtest,f'oceanm_{YR0}_{jday:03d}.nc')
      #print(f'Reading {docn}')
      with xarray.open_dataset(docn) as ds:
        A2d = ds[varnc].isel(zl=0).data.squeeze()
    elif regn == 'ARC':
      pthtest = f'/archive/Dmitry.Dukhovskoy/fre/ARC12/test_ice_relax/ARCphys_expt{expt_nmb:02d}/{YRS}-{mstart:02d}'
      docn = os.path.join(pthtest,f'ocean_daily.nc')
      dref = mtime.datenum([1993,1,1])
      with xarray.open_dataset(docn, decode_times=False) as ds:
        TIME = ds['time'].data + dref
        DTM = np.abs(TIME-dnmb0)
        itime = np.argmin(DTM)
        assert DTM[itime] < 1, f'Given date: {YR}/{MM}/{DD} - Check dates in the arch file {docn}'
        A2d = ds[varnc].isel(time=itime).data.squeeze()
      
    A2d[HH>=0]=np.nan
    A2d = region_lim(A2d, regn)
    if Acell_msk is None:
      Acell_msk = Acell.copy()
      Acell_msk[np.isnan(A2d)] = np.nan

    # Spatial mean:
    Amn = A2d * Acell_msk
    fld_mn = np.nansum(Amn) / np.nansum(Acell_msk)
    
    FLD[irec,iexp] = fld_mn
    print(f'expt {expt_nmb:02d} mean {varnm}={fld_mn:.4f}')
    iexp += 1
    

ECOLR = msisrlx.irlx_tests_colors()


ndays_mnth = 365/12
TM = np.array(TM)
TMd = TM-TM[0]
TMm = 1+TMd/ndays_mnth  # month fractions

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.5, 0.8, 0.4])
hndls = []
lnw0 = 2
for kk in range(nexpts):
  expt_nmb = EXPTS[kk]
  expt_name = msisrlx.irlx_tests_name(expt_nmb, regn)
  fld_mn = FLD[:,kk]
  rlxt = RLXH[kk]
  clr = ECOLR[kk]
  if expt_nmb == 1:
    # Control run
    marker = 'o'
    mksz = 4
  else:
    marker = None
    mksz = None
  
  ln1, = ax1.plot(TMm, fld_mn, 
    linestyle='-', linewidth=lnw0, marker=marker, markersize=mksz, color=clr, label=expt_name)
    #'-', linewidth=2, color=clr, label=f'expt {expt_nmb:02d} {rlxt:03d}hrs')
  hndls.append(ln1)


xtck = [x for x in range(1,13)]
ax1.set_xticks(xtck)
ax1.set_xlim([1,13])
ax1.grid('on')
ax1.set_xlabel('Months')
ax1.set_ylabel(f'{varnm}')
if regn == 'NEPA':
  regn_name = 'NEP10k Arctic'
elif regn == 'NEPB':
  regn_name = 'NEP10k BeringSea'
elif regn == 'ARC':
  regn_name = 'ARC10k inside rlx zone'
sttl = f'IRLX experiments, {varnm} spat.avrg, {regn_name}'
ax1.set_title(sttl)
 
ax2 = plt.axes([0.7, 0.25, 0.25, 0.18])
ax2.legend(handles=hndls, loc='upper right')
ax2.axis('off')



