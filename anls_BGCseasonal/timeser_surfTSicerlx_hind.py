"""
  Show timeseries of monthly surface T/S
  inside the relaxation zone
  to check for the drifts due to ice relaxation

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
import argparse
import pickle

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
import mod_sis2_relax as msisrlx
import mod_rtofs as mrtofs
import mod_hausdorff_distance as mmhd
importlib.reload(mutob)
importlib.reload(msisrlx)

parser = argparse.ArgumentParser()
parser.add_argument("--yrs", help="year to start", type=int, required=True)
parser.add_argument("--yre", help="year to end", type=int)
parser.add_argument("--regn", help="NEPA - Arct.part of NEP only, NEPB - Ber.Sea NEP only", \
                     type=str, required=True)
#parser.add_argument("--varnm", help="temp or salin", type=str, required=True)
args = parser.parse_args()

# Test runs were performed for only 1 year
regn = args.regn if args.regn else None
YRS = args.yrs if args.yrs else None
YRE = args.yre if args.yre else YRS
#varnm = args.varnm if args.varnm else None

mstart = 1  # current test runs all started on Jan 1, 2001
expt_nmb = 2 # hindcast 2
# relax hours:
RLXH = [24]


fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

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

jdm, idm = HH.shape

# Bering Sea - Chukchi Sea masks:
BMsk, AMsk = msisrlx.mask_NEP10k_BerArc(HH,hlat)
#JB,IB = np.where(BMsk==1)

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
  return A2d 


DX, DY = mmom6.dx_dy(hlon, hlat)
Acell  = DX*DY*1.e-6  # km2

hcst_time = 3 # f/csat time interval, months
hcst_interv = np.array([x for x in range(1,12+hcst_time,hcst_time)], dtype=int)


tvarnc = 'tos'
svarnc = 'sos'

nyrs = YRE-YRS+1
nrecs = nyrs*12
Acell_msk = None
TM = []
TFLD = np.zeros((nrecs))
SFLD = np.zeros((nrecs))
irec = -1
for YR in range(YRS,YRE+1):
  for MM in range(1,13):
    DD = 15
    dnmb0 = mtime.datenum([YR,MM,DD])
    print(f"Processing {YR}/{MM}/{DD} ...")
    irec += 1

    TM.append(dnmb0)

    # Find init date for given month, assuming hcst_time (n months) f/cast interval
    kint = np.searchsorted(hcst_interv, MM, side='right') - 1
    assert(hcst_interv[kint] <= MM < hcst_interv[kint+1]), f'Wrong time bin {kint} for {MMA}'
    MINIT = hcst_interv[kint]
    imo = MM-MINIT      # current month in the archive output

    pthhnd = f'/archive/Dmitry.Dukhovskoy/fre/NEP/hindcast_bgc/NEPbgc_nudged_hindcast02/history/{YR}{MINIT:02d}01'
    docn = os.path.join(pthhnd,f'ocean_month.nc')
    #print(f'Reading {docn}')
    with xarray.open_dataset(docn) as ds:
      T2d = ds[tvarnc].isel(time=imo).data.squeeze()
      S2d = ds[svarnc].isel(time=imo).data.squeeze()
    
    T2d[HH>=0]=np.nan
    T2d = region_lim(T2d, regn)
    S2d[HH>=0]=np.nan
    S2d = region_lim(S2d, regn)
    if Acell_msk is None:
      Acell_msk = Acell.copy()
      Acell_msk[np.isnan(T2d)] = np.nan

    # Spatial mean:
    Tmn = T2d * Acell_msk
    t_mn = np.nansum(Tmn) / np.nansum(Acell_msk)
    Smn = S2d * Acell_msk
    s_mn = np.nansum(Smn) / np.nansum(Acell_msk)
    
    TFLD[irec] = t_mn
    SFLD[irec] = s_mn
    print(f'NEP10k hcast, {regn}  mean SST={t_mn:.4f} SSS={s_mn:.4f}')
    

CLR = [[0.,0.3,1],
       [0.9,0.5,0],
       [0.,0.9,0.7],
       [1.,0.9,0],
       [0.8,0.,0.5],
       [0.7, 1, 0.2]]

TM = np.array(TM)
DV = mtime.datevec1D(TM, fHR=False)
DV = np.array(DV).transpose()
time_yrs = DV[:,0]+(DV[:,1]-1)/12

rlxt = 24
clr = CLR[0]
if regn == 'NEPA':
  regn_name = 'NEP10k Arctic Chukchi'
elif regn == 'NEPB':
  regn_name = 'NEP10k BeringSea'
expt_name = 'NEPbgc_nudged_hindcast02'
xtck = [x for x in range(YRS,YRE+1)]

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.55, 0.8, 0.4])
ax1.plot(time_yrs,SFLD, '-', linewidth=2, color=clr)

ax1.set_xticks(xtck)
ax1.set_xlim([YRS,YRE+1])
ax1.grid('on')
sttl = f'{expt_name} irlx={rlxt} hrs, SSS spat.avrg, {regn_name}'
ax1.set_title(sttl)

# Plot T
ax2 = plt.axes([0.1, 0.08, 0.8, 0.4])
ax2.plot(time_yrs,TFLD, '-', linewidth=2, color=[0.,0.9,0.3])

ax2.set_xticks(xtck)
ax2.set_xlim([YRS,YRE+1])
ax2.grid('on')
sttl = f'{expt_name} irlx={rlxt} hrs, SST spat.avrg, {regn_name}'
ax2.set_title(sttl)

btx = 'timeser_surfTSicerlx_hind.py'
bottom_text(btx, pos=[0.05,0.02])



