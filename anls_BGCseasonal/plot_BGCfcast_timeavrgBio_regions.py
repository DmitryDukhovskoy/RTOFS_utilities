"""
    Plot 2D mean BGC fields 
  BGC seasonal forecasts

  WOA23, figures:
  https://www.ncei.noaa.gov/access/world-ocean-atlas-2023f/

Dissolved Oxygen (µmol/kg)
Percent Oxygen Saturation (%)
Apparent Oxygen Utilization (µmol/kg)
Silicate (µmol/kg)
Phosphate (µmol/kg)
Nitrate (µmol/kg)

  stereographic projection
  in different regions: # CalCur - Calif Current region, GulfAlaska - Alaska region, BeringChuk - Bering
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
import mod_colormaps as mclrmps
importlib.reload(mutob)
importlib.reload(manseas)

# Saved on coarser z:
# z_l = 2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 55, 65, 75, 
#    85, 95, 105, 115, 125, 135, 145, 162.5, 187.5, 212.5, 237.5, 262.5, 
#    287.5, 325, 375, 425, 475, 550, 650, 750, 850, 950, 1050, 1150, 1250, 
#    1350, 1450, 1625, 1875, 2125, 2375, 2750, 3250, 3750, 4250, 4750, 5250, 
#    5750, 6250

parser = argparse.ArgumentParser()
parser.add_argument("--varnm", help="o2 no3 po4 sio4", type=str, required=True)
parser.add_argument("--regnm", help="CalCur BeringChuk GulfAlaska", type=str)
parser.add_argument("--yrs", help="Year to start averaging", type=int, required=True)
parser.add_argument("--yre", help="Year to end averaging, defualt = yrs", type=int)
parser.add_argument("--mms", help="Calendar month, start of avrg", type=int, required=True)
parser.add_argument("--mme", help="Calendar month, end of avrg, default=yrs", type=int)
parser.add_argument("--minit", help="Initialization month 1,4,7,10", type=int, required=True)
parser.add_argument("--ensS", help="Start ens avergaing: Ensemb numb, 1,...,10", type=int, required=True)
parser.add_argument("--ensE", help="End ens avergaing, default=ensS no averaging", type=int)
parser.add_argument("--zz", help="approximate depth to plot: 0, ..., 5000 m", type=float, required=True)
args = parser.parse_args()

varnm  = args.varnm if args.varnm else None
regn_name = args.regnm if args.regnm else 'CalCur'
YRS    = args.yrs if args.yrs else None
YRE    = args.yre if args.yre else YRS
MMS    = args.mms if args.mms else None
MME    = args.mme if args.mme else MMS
MMI    = args.minit if args.minit else None
ensS   = args.ensS if args.ensS else None
ensE   = args.ensE if args.ensE else ensS
zz_plt = args.zz if args.zz is not None else None # ocean depth to plot

zz_plt = -abs(zz_plt)
varnc = varnm
if varnm == 'sio4':
  varnc='si'

YAVRG = [x for x in range(YRS,YRE+1)]
MAVRG = [x for x in range(MMS,MME+1)]
ENSMB = [x for x in range(ensS,ensE+1)]

expt_name = 'NEPbgc_fcst_dailyOB01'
run_info = f'{expt_name} avrg {varnm} MI={MMI} {min(YAVRG)}-{max(YAVRG)} '+\
           f'Mo: {min(MAVRG)}-{max(MAVRG)} e{ensS:02d}-e{ensE:02d}'

print(f'Plotting {varnm} {expt_name} e-{ensS:02d}-e{ensE:02d} ')
print(f'{run_info}')


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

# Read vertical layers from oceanm archive file:
ptharch = f'/archive/Dmitry.Dukhovskoy/fre/NEP/forecast_bgc/{expt_name}'

def read_mnth_mom6(YRI, MMI, ens_nmb, varnc, iselZ=None, iselT=None):
  """
    3 month runs are assumed for the hindcasts
  """
  pthfull = os.path.join(ptharch,f'{YRI}-{MMI:02d}-e{ens_nmb:02d}/history')
  flnm = 'ocean_cobalt_tracers_month_z.nc'
  dflnm = os.path.join(pthfull,flnm)
  print(f'Reading {dflnm}')
  with xarray.open_dataset(dflnm) as ds_mom:
    if (iselZ is not None) and (iselT is not None):
      AA = ds_mom[varnc].isel(time=iselT,z_l=iselZ).data.squeeze()
    else:
      AA = ds_mom[varnc].data

  return AA

ZM = read_mnth_mom6(YRS, MMI, ensS, 'z_l')
ZM = -abs(ZM)
dZ = np.abs(ZM-zz_plt)
ilr0 = np.argmin(dZ)
lr0  = ilr0+1
zz0 = ZM[ilr0]  # actual depth to be plotted


rmin, rmax, tscntrs, tslabels = manseas.colormap_params(regn_name, varnm, zz0=zz0)

Ncmp = 200
log_scale = False
logstr = ''
match varnm:
  case 'o2':
    cff = 1.e6
    unts = 'mcromol/kg'
    clrmp = mclrmps.colormap_haline2(end_clr=[0.7,0.1,0], start_clr=[0.8,0.8,1])

  case 'po4':
    cff = 1.e6
    unts = 'mcromol/kg'
    #clrmp = mclrmps.colormap_temp_coldhot()
    #clrmp = mclrmps.colormap_speed()
    #clrmp = mclrmps.colormap_ice_thkn()
    clrmp = mclrmps.colormap_conc()

  case 'sio4':
    cff = 975.61          # output is mole/m3 --> micro mole/kg, assuming rho_ocean = 1025 kg/m3
    unts = 'mcromol/kg'
    #clrmp = mclrmps.colormap_temp_coldhot()
    #clrmp = mclrmps.colormap_speed()
    clrmp = mclrmps.colormap_ice_thkn()
    #clrmp = mclrmps.colormap_conc()

  case 'no3':
    cff = 1.e6
    log_scale = True
    log_str = 'log'
    unts = 'mcromol/kg'
    CLRS = [[1, 1, 1],
        [0.6, 0.02, 0.6],
        [0.2, 0.38, 1],
        [0., 0.8, 0.8],
        [0.4, 0.8, 0],
        [1, 1, 0.5],
        [1, 0.8, 0.6],
        [1, 0.6, 0],
        [0.7, 0.1, 0.1]]

    #clrmp = mclrmps.colormap_temp(clr_ramp=[1,1,1])
    clrmp = mclrmps.colormap_posneg_uneven(CLRS)

clrmp.set_bad(color=[0., 0., 0.])


Time = []
irec = 0
# 1995 - 2004
for YR0 in YAVRG:
  for MM0 in MAVRG:
    dnmb0 = mtime.datenum([YR0, MM0, 15])
    YRI = manseas.yr_init_fcst_from_datenum(dnmb0, MMI)  # init year
    if YRI < 1993:
      continue    # cycle, outside the f/cast time period
    elif YRI > 2024:
      continue

    Time.append(dnmb0)

    for ens_nmb in ENSMB:
      MMF = manseas.mofcst_from_mocalend(YRI, MMI, MM0)
      imo = MMF-1      # forecast month index in the archive output
      print(f'Reading ens {ens_nmb:02d}, MMI={MMI} MMF={MMF} imo={imo}')
      AA  = read_mnth_mom6(YRI, MMI, ens_nmb, varnc, iselZ=ilr0, iselT=imo)

      if irec == 0:
        A2d = AA.copy()
      else:
        A2d = A2d + AA

      irec += 1

A2d = A2d/irec * cff
DV = mtime.datevec1D(Time)

# Mask ocean > zmin depth:
#A2d = np.where( (np.isnan(A2d)) & (HH<0), -1.e3, A2d)
if lr0 > 2:
  A2d = np.where(HH>=zz0, np.nan, A2d)


if log_scale:
  JJ0,II0 = np.where(A2d <= 1.e-32)
  if len(JJ0) > 0:
    A2d[JJ0,II0] = np.nan
  lA2d = np.log(A2d)
  if len(JJ0) > 0:
    lA2d[JJ0,II0] = 0.

  A2d = lA2d.copy()


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


btx = 'plot_BGCfcast_timeavrgBio_regions.py'
sttl = f"{run_info} z={zz0:6.1f} m"

# Stereographic projection:
from mpl_toolkits.basemap import Basemap, cm

lon0, lat0, height, width = manseas.stereogr_params_regions(regn_name)

m = Basemap(width=width, height=height, resolution='l',\
            projection='stere', lat_ts=55, lat_0=lat0, lon_0=lon0)

xR, yR = m(hlon, hlat)

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()

ax1 = manseas.plot_stereogr_axis(fig1, m, xR, yR, A2d, clrmp, rmin, rmax, \
                       btx=btx, tscntrs=tscntrs, tslabels=tslabels, sttl=sttl)

# Plot NEP domain:
plt.sca(ax1)
xdom, ydom = m(Xreg, Yreg)
m.plot(xdom, ydom, 'w-')





