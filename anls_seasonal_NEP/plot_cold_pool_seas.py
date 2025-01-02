"""
  Plot bottom T and cold pool in the Bering Sea by seasons:
  bottom water < 2C
  e.g. On the variability of the Bering Sea Cold Pool and implications 
       for the biophysical environment
  2022
 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8979450/

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
importlib.reload(manseas)


# Initial date
# Look at ens run #1 - the only ens. that has 5-day av. output fields
expt     = 'seasonal_daily'  # seasonal forecasts with dailyOB from SPEAR
varnm    = 'salin'  # temp (potential) / salin
#dnmbS    = mtime.datenum([2015,1,1])
# Averaging time period:
MMS   = 1    # f/cast init. month in each year, can be changed to months: 1, 4, 7, 10
YAVRG = [x for x in range(2005,2015)]
#MAVRG = [1,2,3]  # months to average: Winter  JFM, Summer: JAS
#MAVRG = [4,5,6]  # months to average: Spring, AMJ
#MAVRG = [7,8,9]  # months to average: Winter  JFM, Summer: JAS
MAVRG = [10,11,12]  # months to average: Fall

nensR    = 1
expt_nmb = 2   # 2 - seas f/casts with dailyOB


expt_name = f'NEPphys_frcst_dailyOB-expt{expt_nmb:02d}'
run_info = f'{expt_name} init MM={MMS} e{nensR:02d}, conservT bottom: {min(YAVRG)}-{max(YAVRG)} Mo: {min(MAVRG)}-{max(MAVRG)}'

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

# Get indices of the polygon:
II = pthseas['ANLS_NEP']['poly_BerSea']['II']
JJ = pthseas['ANLS_NEP']['poly_BerSea']['JJ']
jdm, idm = HH.shape

#DX, DY = mmom6.dx_dy(hlon, hlat)
#Acell  = DX*DY
X, Y   = np.meshgrid(np.arange(idm), np.arange(jdm))
MS, _, _ = mmisc.inpolygon_v2(X, Y, II, JJ)  # 
JBS, IBS = np.where( (MS == 1) & (HH >= -250) & (HH < 0) ) #exclude deeep regions
MSKBS  = np.zeros((jdm,idm))
MSKBS[JBS,IBS] = 1


# Save contour of the Bering Sea shelf for WOA23:
pthanls = pthseas['MOM6_NEP'][expt]['pthanls'].format(expt_nmb=expt_nmb) 
dfnm_rbnd = os.path.join(pthanls,f'BeringShelf_boundary_{jdm}x{idm}.pkl')
dfnm_vrtx = os.path.join(pthanls,f'BeringShelf_region_vrtx.pkl')
msk_save = not os.path.isfile(dfnm_rbnd)
if msk_save:
  import pickle
  print('Deriving reigonal boundary for other analyses')
  yC = int(np.mean(np.where(MSKBS==1)[0]))
  xC = int(np.mean(np.where(MSKBS==1)[1]))
  RBOUND = manseas.derive_closed_contour(MSKBS, xC, yC, hcntr=[0.9], indx_int=True)
  RBOUND = manseas.smooth_coastline(RBOUND, npnts=11, indx_int=True)
  Ibnd = RBOUND[:,0]
  Jbnd = RBOUND[:,1] 
  Xbnd = hlon[Jbnd,Ibnd]
  Ybnd = hlat[Jbnd,Ibnd]
  # Save for analysis using WOA23 or other fields:
  print(f'Dumping Bering Shelf bndry --> {dfnm_rbnd}')
  with open(dfnm_rbnd, 'wb') as fid:
    pickle.dump([Xbnd,Ybnd],fid)

  # Save original contour pnts:
  XBS = hlon[JJ,II]
  YBS = hlat[JJ,II]
  print(f'Dumping Bering Shelf bndry --> {dfnm_vrtx}')
  with open(dfnm_vrtx, 'wb') as fid:
    pickle.dump([XBS,YBS],fid)

#mo_fcsts = manseas.yrmo_seasonal_fcst(YRS, MMS)

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
dP = manseas.read_oceanm3D_field(pthfull, floceanm, 'h', notime=False)
dP = np.where(dP < 1.e-3, 0., dP)
ZZ = mmom6.zm2zz(ZM)

Time = []
iyr  = 0
for YRS in (YAVRG):
  dnmbS    = mtime.datenum([YRS,MMS,1])
  pthfcst0 = os.path.join(pthoutp0,f'{YRS}-{MMS:02d}-e{nensR:02d}','history')
  TT, TM = manseas.monthly_mean_from_Ndaily_ocean3D(pthfcst0, YRS, MMS, 'temp', ocnfld, MAVRG=MAVRG)
  SS, TM = manseas.monthly_mean_from_Ndaily_ocean3D(pthfcst0, YRS, MMS, 'salt', ocnfld, MAVRG=MAVRG)

  if iyr == 0:
    T3d = TT.copy()
    S3d = SS.copy()
  else:
    T3d = T3d + TT
    S3d = S3d + SS
  Time = Time + TM

  iyr += 1

T3d = T3d/iyr
S3d = S3d/iyr
DV = mtime.datevec2D(Time)

sys.path.append(PPTHN + '/TEOS_10/gsw')
sys.path.append(PPTHN + '/TEOS_10/gsw/gibbs')
sys.path.append(PPTHN + '/TEOS_10/gsw/utilities')
import mod_swstate as msw
import conversions as gsw
# Compute absolute salinity from practical S:
print('Computing absolute S')
jdm, idm = HH.shape
kdm   = len(ZM)
Z3d   = np.tile(ZM, idm*jdm).reshape((idm,jdm,kdm))
Z3d   = np.transpose(Z3d, (2, 1, 0))
PR    = np.zeros((kdm,jdm,idm))
for kk in range(kdm):
  pr_db, _ = msw.sw_press(Z3d[kk,:,:].squeeze(), hlon)
  PR[kk,:] = pr_db

SA = gsw.SA_from_SP(S3d, PR, hlon, hlat)

# Compute conservative T from potential T
print('Computing conservative T')
CT3d = gsw.CT_from_pt(SA, T3d)

# Derive bottom T:
kdm, jdm, idm = T3d.shape
Tbtm = np.zeros((jdm,idm))*np.nan
dpmin = 1.e-1
for ik in range(1,kdm):
  dpup  = dP[ik-1,:].squeeze()
  dpbtm = dP[ik,:].squeeze()
  tz    = CT3d[ik-1,:]
  if ik < kdm-1:
    Jb, Ib = np.where( (dpup > dpmin) & (dpbtm <= dpmin) )
  else:
# Deep layers include all left:
    Jb, Ib = np.where( dpup > dpmin )
  if len(Jb) == 0: continue
  Tbtm[Jb, Ib] = tz[Jb, Ib]

# Mask outside region:
Tbtm = np.where( (MSKBS==0) & (HH<0) , 1.e3, Tbtm)
# Check, should be empty:
j0,i0 = np.where( (np.isnan(Tbtm)) & (HH < -10) )
if len(j0) > 0:
  print(f'WARNING: {len(j0)} points Bottom T is missing')

CLRS = [[0.6, 0.02, 0.6],
        [0.2, 0.38, 1],
        [0., 0.8, 0.5],
        [0.9, 0.6, 0],
        [1, 1, 1]]

clrmp = mclrmps.colormap_posneg_uneven(CLRS)
clrmp.set_bad(color=[0.6,0.6,0.6])

rmin = -1.8
rmax = 2.
tscntrs = [x/10 for x in range(-20,80,5)]
tslabels = [x/10 for x in range(-20,80,5)]

sinfo = 'Conservative T in the near-bottom layer'

# Stereographic Map projection:
from mpl_toolkits.basemap import Basemap, cm
width  = 2200*1.e3
height = 2200*1.e3
lat0   = 62.
lon0   = -172.

m = Basemap(width=width, height=height, resolution='l',\
            projection='stere', lat_ts=55, lat_0=lat0, lon_0=lon0)

xR, yR = m(hlon, hlat)

btx = 'plot_cold_pool_seas.py'
sttl = run_info

plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()

ax1 = manseas.plot_stereogr_axis(fig1, m, xR, yR, Tbtm, clrmp, rmin, rmax, \
                       btx=btx, tscntrs=tscntrs, tslabels=tslabels, sttl=sttl)

plt.sca(ax1)
ax1.contour(xR,yR, Tbtm, [2], linestyles='solid', colors=[(1., 0.4, 0.9)])

ax1.contour(xR,yR,HH,[x for x in range(-8000,0,500)], linestyles='solid', colors=[(0.9,0.9,0.9)], linewidths=1)
# Show region:
ax1.contour(xR,yR, MSKBS, [0.9], linestyles='solid', colors=[(0.8,0.2,0)])

# Plot NEP domain:
#plt.sca(ax1)
#xdom, ydom = m(Xreg, Yreg)
#m.plot(xdom, ydom, 'w-')



