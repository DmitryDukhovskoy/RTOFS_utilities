"""
  Plot seasonal mean and STD ellipses of mean wind stress vectors
  projected onto along-coast axis
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
#import torch
import sys
import pdb
import netCDF4
import importlib
from netCDF4 import Dataset as ncFile
import timeit
#import pickle
import xarray
#import yaml
from yaml import safe_load
from matplotlib.patches import Polygon

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
import mod_colormaps as mcmp
import mod_mom6 as mom6util
import mod_misc1 as mmsc1
import mod_anlsnep as manlsnep
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_anls_seas as manseas
import mod_utils_ob as mutob
importlib.reload(mutob)
importlib.reload(manseas)

expt     = 'seasonal_daily'  # seasonal forecasts with dailyOB from SPEAR
varnm    = 'tau'  # 
# Averaging time period:
MMS   = 1    # f/cast init. month in each year, can be changed to months: 1, 4, 7, 10
YAVRG = [x for x in range(2011,2021)]
MAVRG = [1,2,3]  # months to average: Winter  JFM, Summer: JAS
#MAVRG = [7,8,9]  # months to average: Winter  JFM, Summer: JAS
regn_name = 'CalCur' # CalCur - Calif Current region, Alaska - Alaska region, BerSea - Bering
                     # Following Stoke et al., 2015
hcntr = -20. # approximate isobath to follow for transect
plt_sct = False  # plot map showing transect line

# Average over a depth range:
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

DX, DY = mhycom.dx_dy(hlon, hlat)

II = pthseas['ANLS_NEP'][regn_name]['II']
JJ = pthseas['ANLS_NEP'][regn_name]['JJ']
xlim1 = min(II)
xlim2 = max(II)
ylim1 = min(JJ)
ylim2 = max(JJ)


# Derive grid indices on the shelf:
#z0 = -800.
#ShIJ, IJl, IJw = manlsnep.find_NWP_shelf(HH, hlon, hlat, DX, DY, z0)
#Ldst = ShIJ.Dst
SH  = manseas.derive_NEPcoastline(HH, hcntr=hcntr, indx_int=True, xFS=293, yFS=47)
#SH  = manseas.derive_NEPcoastline(HH, hcntr=-300, indx_int=True, xFS=293, yFS=47)
SH  = manseas.chop_contour(SH, [293, 47], [255, 351])
SHf = manseas.smooth_coastline(SH, npnts=11, indx_int=True)
# Make sure that all points are in ocean grid points to avoid NaN's:
SHf = manseas.costaline_ocean_points(SHf, HH)

# Save only contour within the region:
Ish = SHf[:,0]
Jsh = SHf[:,1]
#Iout = np.where( (Jsh>ylim2) | (Jsh<ylim1) | (Ish>xlim2) | (Ish<xlim1) )[0]
#Ish = np.delete(Ish, Iout)
#Jsh = np.delete(Jsh, Iout)
Hbtm = HH[Jsh, Ish]
Ldist_sh, Ldx_sh = manseas.distance_segments(DX, DY, Ish, Jsh)
Xsh = hlon[Jsh,Ish]
Ysh = hlat[Jsh,Ish]

CST = manseas.derive_NEPcoastline(HH, indx_int=True, xFS=299, yFS=43)
CST  = manseas.chop_contour(CST, [296, 42], [264, 359])
CSTf = manseas.smooth_coastline(CST, npnts=21)
Icst = CSTf[:,0]
Jcst = CSTf[:,1]

# Compute average U/V along the transect
ocnfld = 'ocean_monthly'
pthoutp0 = pthseas['MOM6_NEP'][expt]['pthoutp'].format(expt_nmb=expt_nmb)
YR=2011
MM=4
#subdir=f'oceanm_{YR}{MM:02d}'
pthfcst0 = os.path.join(pthoutp0,f'{YR}-{MM:02d}-e01','history')
#list_files = manseas.list_oceanice_files(pthfcst0, prefix=ocnfld)
pthfull = pthfcst0
floceanm = 'ocean_month.nc'

# Not-rotated vectors:
TX  = []
TY  = []
# Rotated to along-coastline axis
TXr = []
TYr = []
icc  = 0
for YRS in (YAVRG):
  pthfcst0 = os.path.join(pthoutp0,f'{YRS}-{MMS:02d}-e{nensR:02d}','history')
  dfoutp = os.path.join(pthfcst0, floceanm)
#  print(f'reading {dfoutp}')
  dset   = xarray.open_dataset(dfoutp)

  for mo in MAVRG:
    print(f'Processing {YRS}/{mo}')
    itx  = mo-1
    dmm  = dset['tauuo'].isel(time=itx).data
    dmf  = mmom6.fill_land3d(dmm)
    Taux = mmom6.collocateU2H(dmf, 'symmetr', f_land0 = False)
    dmm  = dset['tauvo'].isel(time=itx).data
    dmf  = mmom6.fill_land3d(dmm)
    Tauy = mmom6.collocateV2H(dmf, 'symmetr', f_land0 = False)

    T1Dx = Taux[Jsh,Ish]
    T1Dy = Tauy[Jsh,Ish]
    # Find components normal & parallel to the local coastaline direction
    Tcstx, Tcsty = manseas.projectU1d_to_coastline(T1Dx, T1Dy, Ish, Jsh, Icst, Jcst)

    TX.append(T1Dx)
    TY.append(T1Dy)
    TXr.append(Tcstx)
    TYr.append(Tcsty)
    icc += 1
  
TX = np.array(TX)
TY = np.array(TY)
TXr = np.array(TXr)
TYr = np.array(TYr)

# Compute mean tau vector along the coastline:
Tcst_mean = np.nanmean(TYr, axis=0)
Tcst_std  = np.std(TYr, axis=0)


# Plot mean along-coastline vector components +/- 1 STD
import mod_draw_vector as mvec

ipnts = TXr.shape[1]
nskip = 20
Iloc = [x for x in range(nskip,ipnts, nskip)]
vcol = [0,0.4,0.9]
dcol = [1.,0.5,0]
cff  = 3500.
dy   = 15    # tick ends
dltY = 50    # stagger vectors for better visibility 
y0   = -dltY
signy = -1

plt.ion()
fig1 = plt.figure(1,figsize=(12,6))
plt.clf()
ax1 = plt.axes([0.05, 0.2, 0.92, 0.72])
icc = 0
for ikk in Iloc:
  x0 = Ldist_sh[ikk]
  if icc%2 == 0:
    signy = -signy
  icc += 1
  y0 = y0+signy*dltY
  #print(f'ikk={ikk} signy={signy} y0={y0}')
  y1 = y0
  tau = Tcst_mean[ikk]
  x1 = x0 + tau*cff
  bm1, bm2 = mvec.arrow_vertices([x0,y0],[x1,y1], cf_ahd=0.3)
  ax1.plot([x0,x1],[y0,y1], '-', color=vcol, linewidth=3)
  ax1.plot(bm1[:,0],bm1[:,1], '-', color=vcol, linewidth=3)
  ax1.plot(bm2[:,0],bm2[:,1], '-', color=vcol, linewidth=3)

  # Add stdev:
  sgm = Tcst_std[ikk]
  xm1 = x1 - sgm*cff
  xm2 = x1 + sgm*cff
  ax1.plot([xm1,xm2],[y0,y1], '-', color=dcol, linewidth=1.5)
  ax1.plot([xm1,xm1],[y0-dy,y0+dy], '-', color=dcol, linewidth=1.5)
  ax1.plot([xm2,xm2],[y0-dy,y0+dy], '-', color=dcol, linewidth=1.5)


# Plot scale vector:
tau = 0.025 # N/m2
x0 = np.mean(Ldist_sh)
x1 = x0 + tau*cff
y0 = 200
y1 = y0
bm1, bm2 = mvec.arrow_vertices([x0,y0],[x1,y1], cf_ahd=0.3)
ax1.plot([x0,x1],[y0,y1], '-', color=vcol, linewidth=3)
ax1.plot(bm1[:,0],bm1[:,1], '-', color=vcol, linewidth=3)
ax1.plot(bm2[:,0],bm2[:,1], '-', color=vcol, linewidth=3)
vtxt = f'{tau:.3f} N/m2'
ax1.text(x1,y0+30, vtxt)

ax1.axis('scaled')
ax1.set_ylim([-350, 350])
ax1.set_xlim([100, 1.05*max(Ldist_sh)])
ax1.set_xticks(np.arange(0,4500,250))
ax1.grid(axis='x')
ax1.set_yticks([])
ax1.invert_xaxis()
ax1.set_xlabel('Along-shore distance, km')

run_info = f'{expt_name} init MM={MMS}, Along-coast tau N/m2, avrg: {min(YAVRG)}-{max(YAVRG)} Mo: {min(MAVRG)}-{max(MAVRG)}'

sttl = f'{run_info} \n {regn_name} hcntr={hcntr}'
ax1.set_title(sttl)

btx = 'windstress_coast_seas.py'
bottom_text(btx, fsz=8, pos=[0.05, 0.2])


if plt_sct:
  fgnmb = 5
  sctnm = f'cntr {hcntr}m'
  manseas.plot_sect_map(HH, hlat, hlon, Xsh, Ysh, Ldist_sh, fgnmb, sctnm, proj='stere', btx=btx)



