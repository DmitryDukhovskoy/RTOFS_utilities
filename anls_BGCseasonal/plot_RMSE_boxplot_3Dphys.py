"""
  Calc and plot boxplot of 
  RMSE for monthly 3D fields (ocean_month_z.nc)
  for the 1st month when the ensembles are still close

  from BGC forecasts
  for physical variables

  Read monthly mean oceanm_YYYY_MM.nc 
  on original grid

  ssh, tos, sos, ssu, ssv

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

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_misc1 as mmisc
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_anls_seas as manseas
import mod_utils_ob as mutob
import mod_sis2_relax as msisrlx
importlib.reload(manseas)


parser = argparse.ArgumentParser()
parser.add_argument("--varnm", help="field to analyze: temp, salin, uvel, vvel", type=str, required=True)
parser.add_argument("--zz", help="Depth to analyze, m, default - set of depths", type=float)
parser.add_argument("--ys", help="Start year, default=1995", type=int)
parser.add_argument("--ye", help="End year, default=ys", type=int)
args = parser.parse_args()

varnm  = args.varnm if args.varnm else None
YRS    = args.ys if args.ys else 1995
YRE    = args.ye if args.ye else YRS
zz_plt = args.zz if args.zz else None
if zz_plt is not None:
  zz_plt = -abs(zz_plt)

expt_name = 'NEPbgc_fcst_dailyOB01'  # forecast run name

# If zz is not specified, plot several depths:
# Plotted in N subplots
if zz_plt is None:
  ZZP = [-10, -50, -100, -150] 
else:
  ZZP = [zz_plt]

nzz = len(ZZP)

if nzz < 4:
  ncol = nzz
elif nzz == 4:
  ncol = 2
else: 
  ncol = 3
  
nrow = nzz // ncol
assert ncol*nrow == nzz, f'Specified columns/rows ({ncol}/{nrow}) do not match N months {nmnths}'


CLR = [[0., 0.4, 0.9],
       [0., 0.8, 0.5],
       [0.9, 0.6, 0],
       [0.8, 0., 0.5],
       [0.5, 0., 0.9],
       [0.5, 0.4, 0]]


ilr0 = None

fyaml = 'bgc_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

pthtopo    = pthseas['MOM6_NEP']['hindcast']['pthgrid']
fgrid      = pthseas['MOM6_NEP']['hindcast']['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"]['hindcast']["ftopo"]
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


match varnm:
  case 'temp':
    varnc = 'potT'
  case 'salin':
    varnc = 'salt'
  case 'uvel':
    varnc = 'u'
  case 'vvel':
    varnc = 'v'


def find_depth_indx(ZM, zz_plt, ds_ocean):
  ZM = ds_ocean['zl'].data
  ZM = -abs(ZM)
  dZ = np.abs(ZM-zz_plt)
  ilr0 = np.argmin(dZ)
  lr0  = ilr0+1
  zz0 = ZM[ilr0]

  return ilr0, lr0, zz0

INIT_MONTHS = [1,4,7,10]  # init. months
ninit = len(INIT_MONTHS)
MCAL = []
YCAL = []
for yr in range(YRS,YRE+1):
  for mi in INIT_MONTHS:
    if yr == 1993 and mi == 1:
      continue
    MCAL.append(mi)
    YCAL.append(yr)

len_cal = len(MCAL)  
Nens = 10  # total N of ensemble runs
ens_ref = 1 # reference ensemble number wrt RMSE is being computed

itot = -1  # total rec counter
TM  = [] 
RMSE = np.zeros((len_cal, Nens, nzz))*np.nan  # # of initializations per year, N of depths, N of records/file
ZZ0 = []
jinit = -1
for ii in range(len(MCAL)):
  MINIT = MCAL[ii]
  YINIT = YCAL[ii] 
  jinit += 1

  # Find init date for given month, assuming hcst_time (n months) f/cast interval
  #kint = np.searchsorted(hcst_interv, MM, side='right') - 1
  #assert(hcst_interv[kint] <= MM < hcst_interv[kint+1]), f'Wrong time bin {kint} for {MMA}'
  #MINIT = hcst_interv[kint]
  #imo = MM-MINIT      # current month in the archive output

  # Compute RMSE wrt to ensmb1, so keep ensmb 1 as a reference:
  for ens_nmb in range(1,Nens+1):
    # when ens_ref = ens_nmb - RMSE should be 0
    pthmain = f'/archive/Dmitry.Dukhovskoy/fre/NEP/forecast_bgc/{expt_name}'
    pthfcst = os.path.join(pthmain,f'{YINIT}-{MINIT:02d}-e{ens_nmb:02d}')
    dnmb_ref = mtime.datenum([1993,1,1])

    # Depths:
    for kk in range(nzz):
      zz_plt = ZZP[kk]

      MMF = 1
      # Get calendar year for the init month/year and f/cast month:
      YY, MM = manseas.mocalend_from_mofcast(YINIT, MINIT, MMF)
      print(f'Calc RMSE {varnm} {YINIT}-{MINIT:02d}-e{ens_nmb:02d}: {YY}/{MM:02d} zz={zz_plt}m')

      # Get ens ref:
      pthfcst_ref = os.path.join(pthmain,f'{YINIT}-{MINIT:02d}-e{ens_ref:02d}','history')
      docn_ref = os.path.join(pthfcst_ref, f'oceanm_{YY}_{MM:02d}.nc')
      with xarray.open_dataset(docn_ref, decode_times=False) as ds_ref:
        ZM = ds_ref['zl'].data
        ilr0, lr0, zz0 = find_depth_indx(ZM, zz_plt, ds_ref)    
   
        A2d_ref = ds_ref[varnc].isel(zl=ilr0).data.squeeze()

      # Read ens. run:
      pthfcst = os.path.join(pthmain,f'{YINIT}-{MINIT:02d}-e{ens_nmb:02d}','history')
      docn = os.path.join(pthfcst, f'oceanm_{YY}_{MM:02d}.nc')
      with xarray.open_dataset(docn, decode_times=False) as ds_run:
        A2d_run = ds_run[varnc].isel(zl=ilr0).data.squeeze()

      # RMSE
      sqerr = (A2d_ref - A2d_run)**2
      nB = np.count_nonzero(~np.isnan(sqerr))
      rmseBP = np.sqrt(np.nansum(sqerr)/nB)

      # Register:
      jens = ens_nmb-1
      RMSE[jinit,jens,kk] = rmseBP

      if jinit == 0:
        ZZ0.append(zz0)
    

tmo = np.array([x for x in range(1,13)])  # forecast months

# Exclude row with 0 RMSE for ens_ref
RMSE0 = RMSE.copy()
RMSE = np.delete(RMSE, ens_ref-1, axis=1)
# Summarize stats by depths:
idm1, idm2, idm3  = RMSE.shape
RMSE = RMSE.reshape(idm1 * idm2, idm3)


# Plot
CLR = [[0., 0.4, 0.9],
       [0., 0.8, 0.5],
       [0.9, 0.6, 0],
       [0.8, 0., 0.5],
       [0.5, 0., 0.9],
       [0.5, 0.4, 0]]

grp_names = []
STAT = np.zeros((nzz,5))
for kk in range(nzz):
  AA = RMSE[:,kk]
  lprc = 25.
  Amed = np.median(AA)
  Alprc = np.percentile(AA, lprc)
  Auprc = np.percentile(AA, (100-lprc))
  Amin  = np.min(AA)
  Amax  = np.max(AA)

  STAT[kk,:] = np.array([Amed, Alprc, Auprc, Amin, Amax])

  zz0 = ZZ0[kk]
  grp_names.append(f'{abs(np.round(zz0)):.0f}') 
 
sttl = f'RMSE all ens wrt e-{ens_ref:02d} {varnm} {YRS} z={zz0:.2f}m'

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
fig1.clf()  
ax1 = plt.axes([0.1, 0.5, 0.6, 0.4])

ax1 = mutil.plot_boxplot2D(ax1, STAT, sttl=sttl, CLRS=CLR,  \
                    xlbl='Depth, m', ylbl='RMSE', grp_names=grp_names, btx=[])



sinfo = f'{expt_name} RMSE {varnm} during 1st month of the forecast\n' + \
        f'all ensemble runs wrt to {ens_ref:02d}\n'
sinfo = sinfo + f'time period: {YRS}-{YRE} init months: 1,4,7,10, pooled'

ax3 = plt.axes([0.1, 0.3, 0.7, 0.14])
ax3.text(0,0, sinfo)
ax3.axis('off')

btx = 'plot_RMSE_boxplot_3Dphys.py'
bottom_text(btx, fsz=8, pos=[0.1, 0.25])


