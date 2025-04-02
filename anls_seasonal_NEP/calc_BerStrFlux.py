"""
  Calculate monthly mean Bering Strait heat / FW fluxes
  Saved indiviual years in separate files

  usage: run calc_BerStrFlux.py --expt=3 --MMI=7 --YRS=1993 --YRE=2008
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
import pickle
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

parser = argparse.ArgumentParser()
parser.add_argument("--expt", help="experiment number: 1, ...", type=int)
parser.add_argument("--MMI", help="init month of f/cast, 1,4,7,10", type=int)
parser.add_argument("--YRS", help="start year of f/casts to derive cold poot stat, 1993, ...", type=int)
parser.add_argument("--YRE", help="end year of f/casts to derive cold pool stat, 1993, ...", type=int)
args = parser.parse_args()

# Default values: 
YRS    = 1993 # year start of the forecast
YRE    = YRS
MMI    = 4
nens   = 1    # ens # for ensemble runs - =1 for 3D ocean fields
expt_nmb = 3  # =2 - seas. f/casts no irelax, =3 - seas. f/casts with ice relax
Tref   = -1.9   # Ref t for computing heat flux, following Woodgate 2018
Sref   = 34.8 # Ref S for FW flux, following Woodgate 2018


if args.expt:
  expt_nmb = args.expt
if args.YRS:
  YRS = args.YRS
if args.YRE:
  YRE = args.YRE
if args.MMI:
  MMI = args.MMI

expt_name = "seasonal_daily"
runname   = f"NEPphys_frcst_dailyOB-expt{expt_nmb:02d}"

if MMI==1 and YRS==1993:
  print(f"first start year for init month {MMI} is 1994, changing {YRS} to 1994")
  YRS = 1994

fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

pthtopo    = pthseas['MOM6_NEP'][expt_name]['pthgrid']
fgrid      = pthseas['MOM6_NEP'][expt_name]['fgrid']
ftopo_mom  = pthseas["MOM6_NEP"][expt_name]["ftopo"]
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

DX, DY = mmom6.dx_dy(hlon, hlat)
#Acell  = DX*DY

# Define Bering Strait:
IIb = [201, 201]
JJb = [680, 667]   # <--- this will result in positive flux to the Arctic Oc. 
#JJb = [667, 680]  # <--- this will result in positive flux to the Ber. Sea

# Define segment lengths, norms, etc:
# positive: norm vector is to the left as follow the section
# negative: to the right
import mod_mom6_valid as mom6vld
IJ     = np.column_stack((IIb, JJb))
SGMT   = mmisc.define_segments(IJ, DX, DY, curve_ornt='positive', check_pole=False)
II     = SGMT.I_indx
JJ     = SGMT.J_indx
nLeg   = SGMT.Leg_number
LegNorm  = SGMT.LegNorm
#II, JJ = mmisc.xsect_indx(Is, Js)
hLsgm1 = SGMT.half_Lsgm1
hLsgm2 = SGMT.half_Lsgm2
Vnrm1  = SGMT.half_norm1
Vnrm2  = SGMT.half_norm2
XX     = hlon[JJ,II]
YY     = hlat[JJ,II]
Hbtm   = HH[JJ,II]
LSgm   = np.zeros((len(II)))  # total segment length = half1 + half2
for ik in range(len(II)):
   LSgm[ik] = hLsgm1[ik] + hLsgm2[ik]
LSgm[0]  = 2*LSgm[0]
LSgm[-1] = 2*LSgm[-1]

II_hf, JJ_hf, XX_hf, YY_hf, Hb_hf, LSgm_hf = \
           mom6vld.segm_half_coord(II, JJ, hLsgm1, hLsgm2, XX, YY, Hbtm)

# Define weights for U and V for each segment:
nsgm  = len(II)
whtU1 = np.zeros((nsgm))
whtV1 = np.zeros((nsgm))
whtU2 = np.zeros((nsgm))
whtV2 = np.zeros((nsgm))
for isgm in range(nsgm):
  nrm1   = Vnrm1[isgm]  # normal for 1st half of the segment
  nrm2   = Vnrm2[isgm]  # --"-- for 2nd half
  ii0    = II[isgm]
  jj0    = JJ[isgm]
  hlsgm1 = hLsgm1[isgm]
  hlsgm2 = hLsgm2[isgm]
  lnrm1  = np.sqrt(nrm1[0]**2 + nrm1[1]**2)
  lnrm2  = np.sqrt(nrm2[0]**2 + nrm2[1]**2)
# No zero-length norms:
  if lnrm1 < 1.e-30 and lnrm2 < 1.e-30:
    raise Exception(f"segm {isgm} has 0 normal vectors")

  whtU1[isgm] = nrm1[0]
  whtU2[isgm] = nrm2[0]
  whtV1[isgm] = nrm1[1]
  whtV2[isgm] = nrm2[1]


# Derive layer thicknesses:
import mod_mom6 as mom6util
ocnfld = 'oceanm'
pthoutp0 = pthseas['MOM6_NEP'][expt_name]['pthoutp'].format(expt_nmb=expt_nmb)
subdir=f'oceanm_{YRS}{MMI:02d}'
pthfcst0 = os.path.join(pthoutp0,f'{YRS}-{MMI:02d}-e01','history')
list_files = manseas.list_oceanice_files(pthfcst0, prefix=ocnfld, subdir=subdir)
pthfull = os.path.join(pthfcst0,subdir)
floceanm = list_files[0]
ZM = manseas.read_oceanm3D_field(pthfull, floceanm, 'zl', notime=False)
ZM = -abs(ZM)
nlrs = len(ZM)

ZZ = mom6util.zm2zz(ZM)
dZ = abs(np.diff(ZZ))

for YRI in range(YRS,YRE+1):
  icc = 0
  VFlx = np.zeros((12))
  FWFlx = np.zeros((12))
  HFlx = np.zeros((12))
  for MMF in range(12):  
  # Forecasts months, calendar month:
    MCal = np.mod(MMI+MMF,12)
    if MCal == 0: 
      MCal = 12

    pthfcst0 = os.path.join(pthoutp0,f'{YRI}-{MMI:02d}-e01','history')
    U2d, TM = manseas.monthly_vsect_mean_from_Ndaily3D(pthfcst0, YRI, MMI, 'u', \
                                    ocnfld, II, JJ,  MAVRG=[MCal], nlrs=nlrs)

    V2d, _ = manseas.monthly_vsect_mean_from_Ndaily3D(pthfcst0, YRI, MMI, 'v', \
                                    ocnfld, II, JJ,  MAVRG=[MCal], nlrs=nlrs)

    T2d, _ = manseas.monthly_vsect_mean_from_Ndaily3D(pthfcst0, YRI, MMI, 'potT', \
                                    ocnfld, II, JJ,  MAVRG=[MCal], nlrs=nlrs)

    S2d, _ = manseas.monthly_vsect_mean_from_Ndaily3D(pthfcst0, YRI, MMI, 'salt', \
                                    ocnfld, II, JJ,  MAVRG=[MCal], nlrs=nlrs)

    # Compute vol, T, S fluxes
    # Use norm direction to define positive flux !!!
    Unrm1 = U2d.copy()*0.
    Unrm2 = U2d.copy()*0.
    nsgm = len(II)
    for isgm in range(nsgm):
      Unrm1[:,isgm] = U2d[:,isgm]*whtU1[isgm] + \
                      V2d[:,isgm]*whtV1[isgm]
      Unrm2[:,isgm] = U2d[:,isgm]*whtU2[isgm] + \
                      V2d[:,isgm]*whtV2[isgm]

      # Arrays for half segments
      # For plotting 2D Unrm section keep original sign of U and V
      UV1 = U2d.copy()*0.
      UV2 = U2d.copy()*0.
      for isgm in range(nsgm):
        UV1[:,isgm] = U2d[:,isgm]*abs(whtU1[isgm]) + \
                      V2d[:,isgm]*abs(whtV1[isgm])
        UV2[:,isgm] = U2d[:,isgm]*abs(whtU2[isgm]) + \
                      V2d[:,isgm]*abs(whtV2[isgm])


    # Vol transport:
    VFlx1 = mom6util.vol_transp_2Dsection(hLsgm1, ZZ, Unrm1)
    VFlx2 = mom6util.vol_transp_2Dsection(hLsgm2, ZZ, Unrm2)
    # Net transport: combine half grid cells:
    VFlx[icc] = np.nansum(VFlx1 + VFlx2)

    # FW flux:
    fw = (Sref-S2d)/Sref
    FWFlx1 = mom6util.vol_transp_2Dsection(hLsgm1, ZZ, (fw*Unrm1))
    FWFlx2 = mom6util.vol_transp_2Dsection(hLsgm2, ZZ, (fw*Unrm2))
    FWFlx[icc] = np.nansum(FWFlx1 + FWFlx2)

    # Heat flux:
    # Negative flux for Bering if T>Tref (loosing heat)
    # assumed: norm is toward AO
    Cp = 4200. # ocean specific heat capacity
    rho_oc = 1025. 
    fheat = -Cp*rho_oc*(T2d-Tref) # note that sign of heat flux depends on norm definition
    HFlx1 = mom6util.vol_transp_2Dsection(hLsgm1, ZZ, (fheat*Unrm1))
    HFlx2 = mom6util.vol_transp_2Dsection(hLsgm2, ZZ, (fheat*Unrm2))
    HFlx[icc] = np.nansum(HFlx1 + HFlx2)

#    U2d = np.expand_dims(U2d, axis=0)
#    V2d = np.expand_dims(V2d, axis=0)
    T2d = np.expand_dims(T2d, axis=0)
    S2d = np.expand_dims(S2d, axis=0)
    if icc == 0:
#      UU = U2d.copy()
#      VV = V2d.copy()
      TT = T2d.copy()
      SS = S2d.copy()
    else:
#      UU = np.append(UU, U2d, axis=0)
#      VV = np.append(VV, V2d, axis=0)
      TT = np.append(TT, T2d, axis=0)
      SS = np.append(SS, S2d, axis=0)

  # Data for plotting 2D sections
  # Project U on the normal vector for the main section line
    UV2d = np.zeros((nlrs,nsgm))
    for isgm in range(nsgm):
  # indices for 1st half segment 
      uu = U2d[:,isgm]
      vv = V2d[:,isgm]
      Snrm = LegNorm[isgm,:]
      UV2d[:,isgm] = uu*Snrm[0] + vv*Snrm[1]


    icc += 1

  pthanls = pthseas['MOM6_NEP'][expt_name]['pthanls'].format(expt_nmb=expt_nmb)
  dflnm = os.path.join(pthanls,f'mnthly_BerSea_Fluxes_expt{expt_nmb:02d}_{YRI}{MMI:02d}.pkl')
  print(f'Dumping Fluxes  --> {dflnm}')
  with open(dflnm, 'wb') as fid:
    pickle.dump([VFlx,FWFlx,HFlx,UV2d,TT,SS,ZM,ZZ,Hbtm,LSgm,XX,YY], fid)






