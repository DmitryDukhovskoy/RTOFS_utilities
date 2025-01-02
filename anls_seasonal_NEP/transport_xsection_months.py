"""
  Calc. tranposrt Sv in the upper layers through the xsection
  Fields are averaged by months
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

from mod_utils_fig import bottom_text
import mod_plot_xsections as mxsct
import mod_time as mtime
import mod_utils as mutil
import mod_read_hycom as mhycom
import mod_misc1 as mmisc
#import mod_valid_utils as mvutil
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_anls_seas as manseas
import mod_utils_ob as mutob
importlib.reload(mutob)

expt     = 'seasonal_daily'  # seasonal forecasts with dailyOB from SPEAR
varnm    = 'v'  # choose U component normal to the section
#dnmbS    = mtime.datenum([2015,1,1])
# Averaging time period:
MMS   = 1    # f/cast init. month in each year, can be changed to months: 1, 4, 7, 10
YAVRG = [x for x in range(1994,2021)]
MWIN = [1,2,3]  # months to average: Winter  JFM, Summer: JAS
MSUM = [7,8,9]  # months to average: Winter  JFM, Summer: JAS
#sctnm  = 'xsct_BerSea' 
sctnm = 'xsct_CalUnderCur_Concep'
Zavrg = -500. # depth for vertical averaging
plt_section = False # show section on the map

nensR    = 1
expt_nmb = 2   # 2 - seas f/casts with dailyOB

expt_name = f'NEPphys_frcst_climOB{expt_nmb:02d}'
run_info  = f'{expt_name} init MM={MMS} e{nensR:02d}, avrg {varnm}: {min(YAVRG)}-{max(YAVRG)} Mo: JFM/JAS'

print(f'Plotting {sctnm} {varnm} {expt_name} ')
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


# Function to print mouse click event coordinates
def onclick(event):
   print([event.xdata, event.ydata])

Is  = pthseas['ANLS_NEP'][sctnm]['II']
Js  = pthseas['ANLS_NEP'][sctnm]['JJ']

IJ     = np.column_stack((Is, Js))
DX, DY = mmom6.dx_dy(hlon,hlat)
SGMT   = mmisc.define_segments(IJ, DX, DY, curve_ornt='positive', check_pole=False)
II     = SGMT.I_indx
JJ     = SGMT.J_indx
nLeg   = SGMT.Leg_number
#II, JJ = mmisc.xsect_indx(Is, Js)
hLsgm1 = SGMT.half_Lsgm1
hLsgm2 = SGMT.half_Lsgm2
XX     = hlon[JJ,II]
YY     = hlat[JJ,II]
Hbtm   = HH[JJ,II]
LSgm   = np.zeros((len(II)))  # total segment length = half1 + half2
for ik in range(len(II)):
   LSgm[ik] = hLsgm1[ik] + hLsgm2[ik]
LSgm[0]  = 2*LSgm[0]
LSgm[-1] = 2*LSgm[-1]


# Distance from the coast:
# Count distances from the coast line, on land = distance<0
Xdist = np.zeros((len(LSgm)))
for ii in range(len(LSgm)-2,-1,-1):
  Xdist[ii] = Xdist[ii+1] + 0.5*(LSgm[ii]+LSgm[ii+1])*1.e-3 # km
Ilnd = int(min(np.where(Hbtm>=0)[0]))
if not Ilnd:
  # no coast:
  Ilnd = 0
else:
  xoffst = 0.5*(Xdist[Ilnd-1]+Xdist[Ilnd])
  Xdist = Xdist - xoffst

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
ZZ   = mmom6.zm2zz(ZM)

npnts = len(LSgm)
#dx_sgm = LSgm.copy()
dx_sgm = np.expand_dims(LSgm, axis=0)
for ik in range(1,nlrs):
  aa = np.expand_dims(LSgm, axis=0)
  dx_sgm = np.append(dx_sgm, aa, axis=0)

DZ = abs(np.diff(ZZ))
DZ = np.expand_dims(DZ, axis=1)
dz_sgm = DZ.copy()
for ik in range(1,npnts):
  dz_sgm = np.append(dz_sgm, DZ, axis=1)

# Truncate bottom or min depth:
for ik in range(npnts):
  hb = Hbtm[ik]
  if Hbtm[ik] < Zavrg:
    hb = Zavrg
  if hb >= 0.:
    dz_sgm[:,ik] = 0.
    continue
  iz = min(np.where(ZZ < hb)[0])
  dz_sgm[iz-1,ik] = abs(ZZ[iz-1] - hb)
  dz_sgm[iz:,ik] = 0.

Time = []
iyr = 0
TRSUM = []
TRWIN = []
for YRS in (YAVRG):
  dnmbS    = mtime.datenum([YRS,MMS,1])
  pthfcst0 = os.path.join(pthoutp0,f'{YRS}-{MMS:02d}-e{nensR:02d}','history')
  AA, _ = manseas.monthly_vsect_mean_from_Ndaily3D(pthfcst0, YRS, MMS, varnm, \
                                    ocnfld, II, JJ,  MAVRG=MSUM, nlrs=nlrs)
  # Depth integrate:
  Udpth_intgr = np.sum(AA*dx_sgm*dz_sgm, axis=0)

  TRSUM.append(Udpth_intgr)

  AA, _ = manseas.monthly_vsect_mean_from_Ndaily3D(pthfcst0, YRS, MMS, varnm, \
                                    ocnfld, II, JJ,  MAVRG=MWIN, nlrs=nlrs)
  # Depth integrate:
  Udpth_intgr = np.sum(AA*dx_sgm*dz_sgm, axis=0)

  TRWIN.append(Udpth_intgr)

  Time.append([YRS])

nyrs  = len(YAVRG)
TRSUM = np.array(TRSUM)
TRWIN = np.array(TRWIN)
#TRNSP = np.reshape(TRNSP,(nyrs,12,npnts))
Time  = np.array(Time)

prct = 25.
TrS_mean = np.nanmean(TRSUM, axis=0)*1.e-6 # Sv
TrS_min  = np.percentile(TRSUM, prct, axis=0)*1.e-6 
TrS_max  = np.percentile(TRSUM, (100-prct), axis=0)*1.e-6

TrW_mean = np.nanmean(TRWIN, axis=0)*1.e-6 # Sv
TrW_min  = np.percentile(TRWIN, prct, axis=0)*1.e-6 
TrW_max  = np.percentile(TRWIN, (100-prct), axis=0)*1.e-6

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.08, 0.3, 0.8, 0.62])
ln1, = ax1.plot(Xdist, TrS_mean, '-', linewidth=2, color=[0.9, 0.5, 0.], label=f'Summer: {min(MSUM)}-{max(MSUM)}')
ax1.plot(Xdist, TrS_min, '-', linewidth=2, color=[1, 0.85, 0.8])
ax1.plot(Xdist, TrS_max, '-', linewidth=2, color=[1, 0.85, 0.8])

ln2, = ax1.plot(Xdist, TrW_mean, '-', linewidth=2, color=[0., 0.4, 0.9], label=f'Winter: {min(MWIN)}-{max(MWIN)}')
ax1.plot(Xdist, TrW_min, '-', linewidth=2, color=[0.8, 0.9, 1])
ax1.plot(Xdist, TrW_max, '-', linewidth=2, color=[0.8, 0.9, 1])

tmax1 = int(1.2*max([np.nanmax(TrS_max),np.nanmax(TrW_max)])*1.e2)*1e-2
tmax2 = int(1.2*max([np.nanmax(abs(TrS_min)),np.nanmax(abs(TrW_min))])*1.e2)*1e-2
tmax  = max([tmax1, tmax2])
yl1 = -tmax
yl2 = tmax

#ax1.axis('scaled')
ax1.set_xlim([0, max(Xdist)])
ax1.set_yticks(np.arange(-1,1,0.1))
ax1.set_ylim([yl1, yl2])
ax1.grid('on')
ax1.set_xlabel('Offshore distance, km')
ax1.invert_xaxis()

run_info = f'{expt_name} init MM={MMS}, avrg: {YAVRG[0]}-{YAVRG[-1]} Mo: {min(MAVRG)}-{max(MAVRG)}, prct={prct:.1f}'

lat_sct = np.mean(hlat[JJ,II])
sttl = f'{run_info} \n Depth-integr Transport, Sv, z0={Zavrg:.0f}, {sctnm} lat={lat_sct:.2f}'
ax1.set_title(sttl)

ax2 = plt.axes([0.1, 0.08, 0.8, 0.14])
lgd = plt.legend(handles=[ln1,ln2], loc='upper right')
ax2.axis('off')

# Compute mean transports:
ii=np.where(TrW_mean>0)[0]
ii1 = ii[0]-1
TrUCw = np.nansum(TrW_mean[ii])    # Undercurrent winter
TrCCw = np.nansum(TrW_mean[:ii1])  # Calif. Current winter

ii=np.where(TrS_mean>0)[0]
ii1 = ii[0]-1
TrUCs = np.nansum(TrS_mean[ii])
TrCCs = np.nansum(TrS_mean[:ii1])  # Calif. Current summer

sinfo = f'Mean transports:\n'
sinfo = sinfo + 'Winter: \n'
sinfo = sinfo + f'  CC = {TrCCw:.2f} Sv\n'
sinfo = sinfo + f'  UC = {TrUCw:.2f} Sv\n'
sinfo = sinfo + 'Summer: \n'
sinfo = sinfo + f'  CC = {TrCCs:.2f} Sv\n'
sinfo = sinfo + f'  UC = {TrUCs:.2f} Sv\n'

ax3 = plt.axes([0.1, 0.08, 0.5, 0.14])
ax3.text(0,0, sinfo)
ax3.axis('off')

btx = 'transport_xsection_months.py'
bottom_text(btx, fsz=8, pos=[0.05, 0.03])


# Bind the button_press_event with the onclick() method
#  fig1.canvas.mpl_connect('button_press_event', onclick)

if plt_section:
  fig2 = plt.figure(2,figsize=(9,8))
  plt.clf()
  ax22 = plt.axes([0.08, 0.2, 0.8, 0.72])
  Lmsk = np.where(HH<0.,1.,0.)

  if min(JJ) == max(JJ):
    xlim1 = min(II)-80
    xlim2 = max(II)+50
    dx = int((xlim2-xlim1)/2)
    dy = 150
    ylim1 = min(JJ)-dy
    ylim2 = max(JJ)+dy

  # Colormaps:
  cmp_Lmsk = mclrmps.colormap_landmask()
  ax22.pcolormesh(Lmsk, cmap=cmp_Lmsk, vmin=0, vmax=0.2)
  ax22.axis('scaled')
  ax22.contour(HH,[x for x in range(-8000,0,1000)], linestyles='solid', \
              linewidths=1, colors=[(0.6, 0.6, 0.6)])
  ax22.contour(hlon,[x for x in range(150,360,5)], linewidths=1, colors=[(0.9, 0.9, 0.9)])
  ax22.contour(hlat,[x for x in range(10,90,5)], linewidths=1, colors=[(0.9, 0.9, 0.9)])
  ax22.plot(II,JJ,'-', linewidth=2, color=[0.9, 0.4,0])
  ax22.set_xlim(xlim1,xlim2) 
  ax22.set_ylim(ylim1,ylim2) 

  sttl2 = f'{sctnm} lat={lat_sct:.2f}'
  ax22.set_title(sttl2)

  bottom_text(btx, fsz=8, pos=[0.05, 0.08])





