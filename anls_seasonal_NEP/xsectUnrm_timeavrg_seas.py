"""
  Plot sections with vertical distribution of U normal to the section
  Fields are averaged over years for specified months 
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
varnm    = 'v'  # potT (potential) / salt
#dnmbS    = mtime.datenum([2015,1,1])
# Averaging time period:
MMS   = 1    # f/cast init. month in each year, can be changed to months: 1, 4, 7, 10
YAVRG = [x for x in range(2011,2021)]
#MAVRG = [1,2,3]  # months to average: Winter  JFM, Summer: JAS
MAVRG = [7,8,9]  # months to average: Winter  JFM, Summer: JAS
#sctnm  = 'xsct_BerSea' 
sctnm = 'xsct_CalUnderCur_Concep'
plt_section = False # show section on the map

# Average over a depth range:
nensR    = 1
expt_nmb = 2   # 2 - seas f/casts with dailyOB

expt_name = f'NEPphys_frcst_climOB{expt_nmb:02d}'
run_info  = f'{expt_name} init MM={MMS} e{nensR:02d}, avrg {varnm}: {min(YAVRG)}-{max(YAVRG)} Mo: {min(MAVRG)}-{max(MAVRG)}'

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

Time = []
iyr = 0
for YRS in (YAVRG):
  dnmbS    = mtime.datenum([YRS,MMS,1])
  pthfcst0 = os.path.join(pthoutp0,f'{YRS}-{MMS:02d}-e{nensR:02d}','history')
  AA, TM = manseas.monthly_vsect_mean_from_Ndaily3D(pthfcst0, YRS, MMS, varnm, \
                                    ocnfld, II, JJ,  MAVRG=MAVRG, nlrs=nlrs)
  if iyr == 0:
    A2d = AA.copy()
  else:
    A2d = A2d + AA

  Time = Time + TM

  iyr += 1

A2d = A2d/iyr

ZZ   = mmom6.zm2zz(ZM)
A2di = mmom6.fill_bottom(A2d, ZZ, Hbtm)

# Plotting
if varnm == 'salin' or varnm == 'salt':
  clrmp = mclrmps.colormap_haline2()
  clrmp.set_bad(color=[0., 0., 0.])
  rmin = 31.6
  rmax = 34.4
#  clrmp.set_under(color=[0.6, 0.6, 0.6])
elif varnm == 'temp' or varnm == 'potT':
  clrmp = mclrpms.colormap_temp(clr_ramp=[0.9,0.8,1])
  clrmp.set_bad(color=[0.,0.,0.])
elif varnm == 'u' or varnm == 'v':
#  cmp_warm = mclrmps.colormap_warm()
#  cmp_cold = mclrmps.colormap_cold()
#  clrmp = mclrmps.colormap_ssh(cpos=cmp_warm, cneg=cmp_cold, nclrs=200)
  clrmp = mclrmps.colormap_uv()
  rmin = -0.1
  rmax = 0.1


plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.08, 0.2, 0.8, 0.72])
im1 = ax1.pcolormesh(Xdist, ZM, A2di, cmap=clrmp, vmin=rmin, vmax=rmax)

# Patch bottom:
verts = [(np.max(Xdist),-8000),*zip(Xdist,Hbtm),(np.min(Xdist),-8000)]
poly = Polygon(verts, facecolor='0.3', edgecolor='0.3', zorder=5)
ax1.add_patch(poly)

#ax1.axis('scaled')
ax1.set_xlim([min(Xdist), max(Xdist)])
ax1.set_yticks(np.arange(-1000,0,100))
ax1.set_ylim([-1000, 0])
ax1.invert_xaxis()

run_info = f'{expt_name} init MM={MMS}, {varnm}, avrg: {min(YAVRG)}-{max(YAVRG)} Mo: {min(MAVRG)}-{max(MAVRG)}'

lat_sct = np.mean(hlat[JJ,II])
sttl = f'{run_info} \n {sctnm} lat={lat_sct:.2f}'
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

btx = 'xsectUnrm_timeavrg_seas.py'
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





