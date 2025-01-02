# Plot U normal to the along-coast transect
# to show upwelling events
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
varnm    = 'u'  # temp (potential) / salin
#dnmbS    = mtime.datenum([2015,1,1])
# Averaging time period:
MMS   = 1    # f/cast init. month in each year, can be changed to months: 1, 4, 7, 10
YAVRG = [x for x in range(2011,2021)]
MAVRG = [1,2,3]  # months to average: Winter  JFM, Summer: JAS
#MAVRG = [7,8,9]  # months to average: Winter  JFM, Summer: JAS
regn_name = 'CalCur' # CalCur - Calif Current region, Alaska - Alaska region, BerSea - Bering
                     # Following Stoke et al., 2015
hcntr = -400. # approximate isobath to follow for transect
plt_sct = True  # plot map showing transect line

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
z0 = -800.
#ShIJ, IJl, IJw = manlsnep.find_NWP_shelf(HH, hlon, hlat, DX, DY, z0)
#Ldst = ShIJ.Dst
SH  = manseas.derive_NEPcoastline(HH, hcntr=hcntr, indx_int=True, xFS=293, yFS=47)
#SH  = manseas.derive_NEPcoastline(HH, hcntr=-300, indx_int=True, xFS=293, yFS=47)
SH  = manseas.chop_contour(SH, [293, 47], [255, 351])
SHf = manseas.smooth_coastline(SH, npnts=11, indx_int=True)


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

# Save only coastline within the region:
#Iout = np.where( (Jcst>ylim2) | (Jcst<ylim1) | (Icst>xlim2) | (Icst<xlim1) )[0]
#Icst = np.delete(Icst, Iout)
#Jcst = np.delete(Jcst, Iout) 
#Icst = Icst.astype(int)
#Jcst = Jcst.astype(int)

# Compute average U/V along the transect
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
nlrs = len(ZM)

Time = []
iyr = 0
for YRS in (YAVRG):
  dnmbS    = mtime.datenum([YRS,MMS,1])
  pthfcst0 = os.path.join(pthoutp0,f'{YRS}-{MMS:02d}-e{nensR:02d}','history')
  UU, TM = manseas.monthly_vsect_mean_from_Ndaily3D(pthfcst0, YRS, MMS, 'u', \
                                    ocnfld, Ish, Jsh,  MAVRG=MAVRG, nlrs=nlrs)

  VV, _  = manseas.monthly_vsect_mean_from_Ndaily3D(pthfcst0, YRS, MMS, 'v', \
                                    ocnfld, Ish, Jsh,  MAVRG=MAVRG, nlrs=nlrs)

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

# Find components normal & parallel to the local coastaline direction
Ucst, Vcst = manseas.projectU_to_coastline(U2d, V2d, Ish, Jsh, Icst, Jcst)

# Plotting
clrmp = mclrmps.colormap_uv()
rmin = -0.15
rmax = 0.15

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.08, 0.2, 0.8, 0.72])
im1 = ax1.pcolormesh(Ldist_sh, ZM, U2d, cmap=clrmp, vmin=rmin, vmax=rmax)

# Patch bottom:
verts = [(np.min(Ldist_sh),-8000),*zip(Ldist_sh,Hbtm),(np.max(Ldist_sh),-8000)]
poly = Polygon(verts, facecolor='0.6', edgecolor='0.6', zorder=5)
ax1.add_patch(poly)

#ax1.axis('scaled')
ax1.set_xlim([min(Ldist_sh), max(Ldist_sh)])
ax1.set_ylim([-600, 0])
ax1.set_yticks(np.arange(-800,0,50))
ax1.invert_xaxis()

tscntrs=[]
if len(tscntrs) > 0:
  CS = ax1.contour(A2d, tscntrs, colors=[cntr_clr], linestyles='solid', linewidths=1)
  if len(tslabels) > 0:
    ax1.clabel(CS, tslabels, inline=1, fontsize=10)

run_info = f'{expt_name} init MM={MMS}, Unormal_coast, avrg: {min(YAVRG)}-{max(YAVRG)} Mo: {min(MAVRG)}-{max(MAVRG)}'

sttl = f'{run_info} \n {regn_name} hcntr={hcntr}'
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


btx = 'offshore_flow_vsect.py'
bottom_text(btx, fsz=8, pos=[0.05, 0.03])


if plt_sct:
  fgnmb = 5
  sctnm = f'cntr {hcntr}m'
  manseas.plot_sect_map(HH, hlat, hlon, Xsh, Ysh, Ldist_sh, fgnmb, sctnm, proj='stere', btx=btx)



