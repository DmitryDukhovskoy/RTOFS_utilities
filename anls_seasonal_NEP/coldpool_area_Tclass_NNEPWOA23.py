"""
  Plot area cold pool on the Bering Sea shelf
  bottom water < 2C for T classes
  Use regional climatology NNEP WOA23 1/10 degree grid

  For plotting: depth-integrate pot. enthalpy and / total depth

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import netCDF4
from netCDF4 import Dataset as ncFile
import importlib
import yaml
from yaml import safe_load
import pickle

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
sys.path.append(PPTHN + '/TEOS_10/gsw')
sys.path.append(PPTHN + '/TEOS_10/gsw/gibbs')
sys.path.append(PPTHN + '/TEOS_10/gsw/utilities')

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_utils as mutil
import mod_read_hycom as mhycom
import mod_colormaps as mcmp
import mod_mom6 as mmom6
import mod_misc1 as mmisc
import mod_anls_seas as manseas
import mod_plot_xsections as mxsct
import matplotlib as mtplt
import mod_regmom as mregmom
import mod_colormaps as mclrmps
import mod_swstate as msw
import conversions as gsw

YR      = 2008 # indicate any year in the decade need to plot

grd=0.25
if grd==0.25:
  cgrd=4
woa='woa23'

woa_seas = {"13": "Jan-Mar",
            "14": "Apr-Jun",
            "15": "Jul-Spt",
            "16": "Oct-Dec",
            "0": "annual"}

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

# Load Bering Shelf boundary, saved in plot_cold_pool_seas.py:
pthanls = pthseas['MOM6_NEP']['seasonal_daily']['pthanls'].format(expt_nmb=2)  
dfnm_rbnd = os.path.join(pthanls,f'BeringShelf_boundary_816x342.pkl')
if not os.path.isfile(dfnm_rbnd):
  raise Exception(f'Region boundary need to be saved in plot_cold_pool_seas.py {dfnm_rbnd}')

with open(dfnm_rbnd, 'rb') as fid:
  RBND = pickle.load(fid)

Xbnd = RBND[0]
Ybnd = RBND[1]

def read_field(furl,varnm):
  print("Reading {1} from {0}".format(furl,varnm))
  nc=ncFile(furl)
# lookup a variable
  dmm0 = nc.variables[varnm][:].data.squeeze()
  dmm = np.copy(dmm0)
  return dmm

def get_lonlat_ZZ(urlT, tfnm, lon1=150., lon2=255., lat1=10., lat2=82.):
  # Get lon/lat for specified region
  furl = os.path.join(urlT,tfnm)
  ZZ  = read_field(furl,'depth')
  ZZ  = -abs(ZZ)
  latW = read_field(furl,'lat')
  lonW = read_field(furl,'lon')
#  lonW = mmisc.shuffle1D_lon180_to0360(lonW0)
  ix1 = np.argmin(np.abs(lonW-lon1))
  ix2 = np.argmin(np.abs(lonW-lon2))+1
  jx1 = np.argmin(np.abs(latW-lat1))
  jx2 = np.argmin(np.abs(latW-lat2))+1
  lonW = lonW[ix1:ix2]
  latW = latW[jx1:jx2]
  jdm  = len(latW)
  idm  = len(lonW)
  LONW = np.zeros((jdm,idm))
  LATW = np.zeros((jdm,idm))
  for ii in range(idm):
    LATW[:,ii]=latW
  for jj in range(jdm):
    LONW[jj,:]=lonW
  
  return [ix1,ix2,jx1,jx2], LONW, LATW, ZZ

def read_TS(urlT, tfnm, var_read):
  # Read T/S from WOA23:
  furl = os.path.join(urlT,tfnm)
#var_read = 't_an'
  A3d = read_field(furl,var_read)
  A3d = np.where(A3d > 1.e10, np.nan, A3d)
# reshaffle to have -180/180 lon inside the domain
#  lonW0 = read_field(furl,'lon')
#  A3d, _ = mmisc.shuffle3D_lon180(A3d, lonW0)

  return A3d

def derive_bottomT(CT3d, dP):
  kdm, jdm, idm = CT3d.shape
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

  return Tbtm

# Cold Pool area by seasons x T classes:
MPLOT  = [x for x in range(1,13,1)]
nmo    = len(MPLOT)
TCLASS = [-2., -1., 0., 1., 2.]
nTC    = len(TCLASS)
CPA    = np.zeros((nmo,nTC-1))   # Cold pool area: 4 seasons x T classes
def extract_regional_clim():
  icc    = 0
  for MM in MPLOT:
    seas, decade, yr1_dec, yr2_dec = manseas.season_decade_woa(YR, MM, month2season=False)
    urlBase = 'https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa/REGCLIM/NNPv2/DATA/'
    urlT    = f"{urlBase}temperature/netcdf/{decade}/0.10/"
    urlS    = f"{urlBase}salinity/netcdf/{decade}/0.10/"
    tfnm    = f"nnp_{decade}_t{seas:02d}_10.nc"
    sfnm    = f"nnp_{decade}_s{seas:02d}_10.nc"

    if icc == 0:
      IJX, LONW, LATW, ZZ = get_lonlat_ZZ(urlT, tfnm)
      ix1,ix2,jx1,jx2 = IJX
      DX, DY = mmom6.dx_dy(LONW, LATW)
      Acell  = DX*DY
      jdm, idm = LONW.shape
      kdm  = len(ZZ)
      Z3d  = np.tile(ZZ, idm*jdm).reshape((idm,jdm,kdm))
      Z3d  = np.transpose(Z3d, (2, 1, 0))

    A3d = read_TS(urlT, tfnm, 't_an')
    T3d = A3d[:,jx1:jx2,ix1:ix2]
    A3d = read_TS(urlS, sfnm, 's_an')
    S3d = A3d[:,jx1:jx2,ix1:ix2]
    kdm, jdm, idm = S3d.shape

    if icc == 0:
      # Land sea mask for WOA
      iz = 0
      LMsk = T3d[iz,:,:].squeeze()
      LMsk = np.where(np.isfinite(LMsk), -10, 1)

      # Derive Regional mask for Bering Shelf:
      Iregn, Jregn = mmisc.find_closest_indx(Xbnd, Ybnd, LONW, LATW)

      X, Y   = np.meshgrid(np.arange(idm), np.arange(jdm))
      MSKBS, _, _ = mmisc.inpolygon_v2(X, Y, Iregn,Jregn)  # 
      Jout, Iout = np.where(MSKBS<1)
      
      # Derive Pressure and layer thicknesses:
      PR   = np.zeros((kdm,jdm,idm))
      for kk in range(kdm):
        pr_db, _ = msw.sw_press(Z3d[kk,:,:].squeeze(), LATW)
        PR[kk,:] = pr_db
      # Estimate layer thickness:
      dP = manseas.derive_dP_WOA(ZZ, T3d)

  # Convert in situ --> potential T at z=0:
    print('Converting T in situ --> T potential')
    Tp3d = mregmom.insitu2pot_3D(T3d, S3d, Z3d, LATW)

    # Compute conservative T:
    print('Computing absolute Salinity')
    SA = gsw.SA_from_SP(S3d, PR, LONW, LATW)
    # Compute conservative T
    CT3d = gsw.CT_from_pt(S3d, Tp3d)

    Tbtm = derive_bottomT(CT3d, dP)
    # Mask outside region:
    Tbtm = np.where( (MSKBS<1.) & (LMsk<0) , 1.e3, Tbtm)
    # Check, should be empty:
    j0,i0 = np.where( (np.isnan(Tbtm)) & (LMsk<0) )
    if len(j0) > 0:
      print(f'WARNING: {len(j0)} points Bottom T is missing')

    for icl in range(1,nTC):
      t1 = TCLASS[icl-1]
      t2 = TCLASS[icl]

      JBS, IBS = np.where( (MSKBS == 1) & (Tbtm <= t2) & (Tbtm > t1))
      if len(JBS) > 0.:
        CP_area = np.sum(Acell[JBS,IBS])*1e-6 # km2
      else:
        CP_area = 0.
      print(f'TCLASS: {t1:.1f} - {t2:.1f} Area={CP_area*1e-5:.1f} x1e5 km2')

      CPA[icc, icl-1] = CP_area

    icc += 1


#  f_save = True
#  if f_save:
  pthanls = pthseas['MOM6_NEP']['seasonal_daily']['pthanls'].format(expt_nmb=2)
  dflout  = os.path.join(pthanls,f'Bering_coldpoolarea_NNEPWOA23_{yr1_dec}-{yr2_dec}.pkl')
  print(f'Saving NNEP Regional Climatology coldpool area --> {dflout}')
  with open(dflout, 'wb') as fid:
    pickle.dump([CPA,TCLASS],fid)

  return()

seas, decade, yr1_dec, yr2_dec = manseas.season_decade_woa(YR, 1, month2season=False)
pthanls = pthseas['MOM6_NEP']['seasonal_daily']['pthanls'].format(expt_nmb=2)
dflout  = os.path.join(pthanls,f'Bering_coldpoolarea_NNEPWOA23_{yr1_dec}-{yr2_dec}.pkl')

if os.path.isfile(dflout):
  print(f'Loading NNEP Regional Climatology coldpool area --> {dflout}')
  with open(dflout, 'rb') as fid:
    CPA, TCLASS = pickle.load(fid)
else:
  print(f'Not found: {dflout}')
  print('Extracting data ...')
  extract_regional_clim() 

 
cff = 1.e-5
CPA = CPA*cff # km2 x 1e-5

CPT   = np.sum(CPA, axis=1)  # total cold pool area

plt.ion()

CLRS = np.array([[0., 0.2, 0.9],
                [0., 0.8, 1],
                [0.7, 0., 1],
                [0.9, 0.4, 0],
                [0.5, 0.3, 0]])

Xtime = [x for x in range(1,13)]
#Xtime = np.arange(1,15,3)


sttl = f"Regional 1/10 Northern NEP WOA23 coldpool area x {cff} km2 ,  monthly climatology {yr1_dec}-{yr2_dec}"


plt.ion()

fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.5, 0.8, 0.4])
LNS = []
for icl in range(nTC-1):
  clr0 = CLRS[icl, :]
  t1 = TCLASS[icl]
  t2 = TCLASS[icl+1]
  tline = f'{t1:.1f} < t < {t2:.1f}'
  ln1, = ax1.plot(Xtime,CPA[:,icl], linewidth=2, color=clr0, label=tline)
  LNS.append(ln1)
#
# Total area for all classes:
clr0 = [0,0,0]
ln1, = ax1.plot(Xtime, CPT, linewidth=2, color=clr0, label='Total')
LNS.append(ln1)

ax1.set_xticks([x for x in range(1,13)])
ax1.set_yticks([x for x in range(0,10)])
ax1.grid(True)
ax1.set_xlim([0.9, 12.1])
ax1.set_ylim([-0.2, 9])
#ax1.set_xticklabels(tck_lbls)
ax1.set_ylabel(f'km2 x {cff}')
ax1.set_xlabel(f'Months')
ax1.set_title(sttl)

ax3 = plt.axes([0.1, 0.2, 0.6, 0.22])
lgd = plt.legend(handles=LNS, loc='upper right')
ax3.axis('off')

btx = 'coldpool_area_Tclass_NNEPWOA23.py'

bottom_text(btx, pos=[0.1, 0.1])

