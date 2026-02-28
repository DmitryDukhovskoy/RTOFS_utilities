"""
  Fix ice concentration, thickness, snow depth in CICE6 restart file
  at the grid cells where the variables exceed threshold values

  The fields are changed to the specified / default values
  adjusting enthalpy  


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
import logging

PPTHN = None
if 'PPTHN' not in locals() or PPTHN is None:
  cwd = os.getcwd()
  parts = cwd.split(os.sep)
  if 'python' in parts:
    idx = parts.index('python')
    PPTHN = os.sep + os.path.join(*parts[:idx + 1])
  else:
    raise RuntimeError("Directory 'python' not found in current working directory path.")

sys.path.extend([
    os.path.join(PPTHN, 'MyPython', 'hycom_utils'),
    os.path.join(PPTHN, 'MyPython', 'draw_map'),
    os.path.join(PPTHN, 'MyPython'),
    os.path.join(PPTHN, 'MyPython', 'mom6_utils')
])

import mod_time as mtime
import mod_cice6_utils as mc6util

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s: %(message)s'
)

logger = logging.getLogger(__name__)

logger.info("Job started")


def find_adj_icepnts(dlti, aice, i0, j0, puny):
  jdm, idm = aice.shape
  i1 = max(i0 - dlti, 0)
  i2 = min(i0 + dlti, idm - 1)
  j1 = max(j0 - dlti, 0)
  j2 = min(j0 + dlti, jdm - 1)

  A = aice[j1:j2+1, i1:i2+1]
  jmm, imm = np.where(A > puny)

  # No grid points with ice found:
  if jmm.size == 0:
    return np.array([], dtype=int), np.array([], dtype=int)

  jice = jmm + j1
  iice = imm + i1

  return iice, jice

def find_adj_ocnpnts(dlti, aice, i0, j0, puny):
  jdm, idm = aice.shape
  i1 = max(i0 - dlti, 0)
  i2 = min(i0 + dlti, idm - 1)
  j1 = max(j0 - dlti, 0)
  j2 = min(j0 + dlti, jdm - 1)

  A = aice[j1:j2+1, i1:i2+1]
  jmm, imm = np.where(A <= puny)

  # No grid points with no ice found:
  if jmm.size == 0:
    return np.array([], dtype=int), np.array([], dtype=int)

  jocn = jmm + j1
  iocn = imm + i1

  return iocn, jocn

def check_print(aice, aicen_new, vicen_new, vsnon_new, ncat):
  # Check hice(n) as it is caclulated in icepack_therm_vertical.F90
  # hice(n) = vice(n) / aice(n) 
  print(' =========  ICE  =========')
  for k in range(1,ncat+1):
    aice_n = aicen_new[k-1,:].squeeze()
    vice_n = vicen_new[k-1,:].squeeze()
    hice_n = np.divide(vice_n, aice_n, out=np.zeros_like(aice), where=aice_n != 0)
    jmin, imin = np.unravel_index(hice_n.argmin(), hice_n.shape)
    jmax, imax = np.unravel_index(hice_n.argmax(), hice_n.shape)
    print(f"Cat {k}, j={jmin}, i={imin}, min hice(n): {np.nanmin(hice_n)}, "+\
          f"aice(n): {aice_n[jmin,imin]}, vice(n): {vice_n[jmin,imin]}")
    print(f"         j={jmax}, i={imax}, max hice(n): {np.nanmax(hice_n)}, "+\
          f"aice(n): {aice_n[jmax,imax]}, vice(n): {vice_n[jmax,imax]}")

  print(' =========  SNOW =========')
  for k in range(1,ncat+1):
    aice_n = aicen_new[k-1,:].squeeze()
    vsno_n = vsnon_new[k-1,:].squeeze()
    hsno_n = np.divide(vsno_n, aice_n, out=np.zeros_like(aice), where=aice_n != 0)
    jmin, imin = np.unravel_index(hsno_n.argmin(), hsno_n.shape)
    jmax, imax = np.unravel_index(hsno_n.argmax(), hsno_n.shape)
    print(f"Cat {k}, j={jmin}, i={imin}, min hsnow(n): {np.nanmin(hsno_n)}, "+\
          f"aice(n): {aice_n[jmin,imin]}, vsno(n): {vsno_n[jmin,imin]}")
    print(f"         j={jmax}, i={imax}, max hsnow(n): {np.nanmax(hsno_n)}, "+\
          f"aice(n): {aice_n[jmax,imax]}, vsno(n): {vsno_n[jmax,imax]}")

  print(' ========= ICE THICKN m3/m2_ice =========')
  vice = np.sum(vicen_new, axis=0).squeeze()
  aice = np.sum(aicen_new, axis=0).squeeze()
  vsnow = np.sum(vsnon_new, axis=0).squeeze()
  vice_m2ice = np.divide(vice, aice, out=np.zeros_like(aice), where=aice != 0)
  vsnow_m2ice = np.divide(vsnow, aice, out=np.zeros_like(aice), where=aice != 0)
  hice_max = np.nanmax(vice_m2ice)
  hsnow_max = np.nanmax(vsnow_m2ice)
  print(f"Max hice = {hice_max:.2f} m3/m2_ice, max hsnow = {hsnow_max:.2f} m3/m2_ice\n")

  return

def fix_hsnow_ithkn(ds_out, vice_m2ice_max, vice_m2ice_fix, vsnow_m2ice_max, vsnow_m2ice_fix,
                    pthrest, flrst_in, logger):
  # CICE parameters:
  puny      = 1.e-11
  c0        = 0.0
  c1        = 1.0
  c2        = 2.0
  p5        = 0.5
  Lsub      = 2.835e6    # latent heat sublimation fw (J/kg)
  Lvap      = 2.501e6    # latent heat vaporization fw (J/kg)
  Lfresh    = Lsub - Lvap # latent heat of melting of fresh ice (J/kg)
  cp_ice    = 2106.       # specific heat of fresh ice (J/ kg/K)
  rhos      = 330.        # density of snow (kg/m3)
  hs_min    = 1.e-4       # min snow thickness for computing Tsno (m)
  nsal      = 0.407
  msal      = 0.573
  min_salin = 0.1      # threshold for brine pocket treatment
  saltmax   = 3.2        # max S at ice base
  hg        = 1.e20    # bad values, land mask, etc.
  nslyr     = 1   # snow layers
  Tmin      = -100.   # minimum snow T

  dflrst_in = os.path.join(pthrest, flrst_in)
  if not os.path.isfile(dflrst_in):
    raise FileNotFoundError(f"Input restart not found: {dflrst_in}")

  if ds_out is None:
    print(f"Reading restart: {dflrst_in}")
    ds_in = xarray.open_dataset(dflrst_in)
    ds_out = ds_in.copy(deep=True)
    ds_in.close()

  # Insert aice and distribute/reduce proportionally over ice cats.:
  assert nslyr == 1, f"Code needs to be modified for nslyr>1, nslyr={nslyr}"
  aicen = ds_out['aicen'].data  # partial area by cats
  vicen = ds_out['vicen'].data  # ice vol per m2 of grid cell by cats
  vsnon = ds_out['vsnon'].data  # snow vol per m2 of ice area
  qsnon = ds_out['qsno001'].data  # snow enthalpy by cats for 1 snow layer
  Tsfcn = ds_out['Tsfcn'].data  # surface T in each cat.
  ncat, jdim, idim = vsnon.shape

  assert aicen.shape == vicen.shape == vsnon.shape

  sice = {}
  qice = {}
  nilrs = 7   # ice layers
  for i in range(1, nilrs+1):
    varnum = f"{i:03d}" 
    sice[varnum] = ds_out[f"sice{varnum}"].data
    qice[varnum] = ds_out[f"qice{varnum}"].data

  # Aggregated ice partial area:
  aice = np.sum(aicen, axis=0).squeeze()

  # mean hice over ice = m3/m2_ice area:
  vice = np.sum(vicen, axis=0).squeeze()
  vice_m2ice = np.divide(vice, aice, out=np.zeros_like(aice), where=aice != 0)

  # mean hsnow over ice = m3/m2_ice area:
  vsnow = np.sum(vsnon, axis=0).squeeze()
  vsnow_m2ice = np.divide(vsnow, aice, out=np.zeros_like(aice), where=aice != 0)

  hice_bad = vice_m2ice > vice_m2ice_max
  hsnow_bad = vsnow_m2ice > vsnow_m2ice_max

  ni_bad = np.count_nonzero(hice_bad)
  ns_bad = np.count_nonzero(hsnow_bad)

  logger.info(f"Found {ni_bad} bad ice points and {ns_bad} bad snow points")
  logger.info(f"Max hice={np.nanmax(vice_m2ice):.2f} m3/m2_ice, hsnow={np.nanmax(vsnow_m2ice):.2f} m3/m2_ice")
  logger.info(f"max hice={vice_m2ice_max:.2f}, max hsnow={vsnow_m2ice_max:.2f}, fix hice={vice_m2ice_fix:.2f}, fix hsnow={vsnow_m2ice_fix:.2f}")

  if ni_bad == 0 and ns_bad == 0:
    logger.info(f"No correction needed, max ice={np.nanmax(vice_m2ice):.1f}m3/m2ice, max snow={np.nanmax(vsnon):.1f}m3/m2ice")
    return ds_out  

  # Select ice points to fix:
  Jins, Iins = np.where(hice_bad | hsnow_bad)
  npnts = Jins.size

  # Ice thickness cats
  # for kitd = 1 - linear remapping
  # and kcatbound = 0 , the lower cat. thickness values:
  # 0.00, 0.64, 1.39, 2.47, 4.57
  hicat = np.array([0., 0.64, 1.39, 2.47, 4.57, 50.])
  dhi_min = 0.01  # min diff between cat ice thicknesses from 2 adjacent cats
  hcat_indx = np.arange(1,ncat+1)

  # Note qsnon, qice < 0 !
  vsnon_new = vsnon.astype(ds_out['vsnon'].dtype).copy()
  qsnon_new = qsnon.astype(ds_out['qsno001'].dtype).copy()
  aicen_new = aicen.astype(ds_out['aicen'].dtype).copy()
  vicen_new = vicen.astype(ds_out['vicen'].dtype).copy()
  Tsfcn_new = Tsfcn.astype(ds_out['Tsfcn'].dtype).copy()
  qicen_new = {}
  sicen_new = {}
  for i in range(1, nilrs+1):
    varnum = f"{i:03d}"
    qicen_new[varnum] = qice[varnum].astype(ds_out[f"qice{varnum}"].dtype).copy()
    sicen_new[varnum] = sice[varnum].astype(ds_out[f"sice{varnum}"].dtype).copy()

  Tfrz = -1.86243522  # ocean freez. T
  Tsfc_max = -0.1  # max surf temp
  dvol_sum = 0.
  dlti = 2         # N of +/- i , j indices to search for ice pnts around i0,j0
  vitot_min = 0.05    # min total ice vol m3/m2_grid, when ai_old = 0 --> ai_new > 0
  hitot_min = np.min([0.1, vice_m2ice_fix])     # mean ice thickness over ice area: = sum(hice(n)*aice(n)) / sum(aice(n)) 
  hsnow_min = np.min([0.01, vsnow_m2ice_fix])   # min snow thickness for noice --> ice case, this is m3/m2_ice 

  for ipp in range(npnts):
    # Points with bad ice values OR bad snow values
    j0 = Jins[ipp]
    i0 = Iins[ipp]

    #if hi_bad[j0,i0]: 
    # Note hsnow = vsn / aice for aice > 0
    # for cat n: vsn(n) = hsnow(n) * aice(n) 
    ai_old  = aice[j0,i0]       # aggreageted ice partial area 
    ain_old = aicen[:,j0,i0]    # partial areas by cats
    vsn_old = vsnon[:,j0,i0]    # snow volume per unit grid-cell area m2
    vin_old = vicen[:,j0,i0]    # ice volume per unit grid-cell area m2 
    tsfcn_old = Tsfcn[:,j0,i0]  # surf T

    # No change in ice conc:
    ai_new = ai_old
    if ai_new <= puny:
      ai_new = 0.
   
    if ai_new > 1.:
      ai_new = 1.

    # iconc does not change, only hice and hsnow:
    if ai_new < puny:
      ain_new = ain_old * 0.
    else:
      if ai_old < puny:
        # all new iconc in cat 1:
        ain_new = ain_old * 0.
        ain_new[0] = ai_new
      else:
        ain_new = ain_old.copy()

    # Update surface T and ice vol / ice thickness
    tsf_new = None
    iice = jice = iocn = jocn = None
    vin_new = None
    aice_case = None
    if ai_old <= puny and ai_new > puny:
      # Case: no ice --> ice, created ice in the grid cell
      aice_case = "noice2ice"

      # Update Tsfcn
      # Find N closest ice points:
      iice, jice = find_adj_icepnts(dlti, aice, i0, j0, puny) 
      if len(iice) == 0:
        # no ice pnt adjacent to j0,i0:
        tsf_new = Tsfcn[:,j0,i0]*0.0 + Tsfc_max
        vin_new = hitot_min * ain_new
      else: 
        # where ice - <= Tmax, where no ice = Tfrz
        Tsf_adj = Tsfcn[:,jice,iice]
        tsf_new  = np.nanmean(Tsf_adj, axis=1)
        tsf_new  = np.where(tsf_new > Tsfc_max, Tsfc_max, tsf_new)
        tsf_new  = np.where(ain_new < puny, Tfrz, tsf_new) 

        vin_new = hitot_min * ain_new

    elif ai_old > puny and ai_new > puny:
      # Case: ice --> updated ice conc
      aice_case = "ice2ice"
      tsf_new = np.where(tsfcn_old > Tsfc_max, Tsfc_max, tsfcn_old)
      tsf_new  = np.where(ain_new < puny, Tfrz, tsf_new)

      hitot_new = np.sum(vin_old)/ai_new
      vin_new   = np.min([vice_m2ice_max, hitot_new]) * ain_new  # m3/m2_cell or grid-cell mean ice thickness, m

    elif ai_old > puny and ai_new < puny:
      # Case: ice --> no ice
      aice_case = "ice2noice"
      # Copy tsfc from adj grid cells
      iocn, jocn = find_adj_ocnpnts(dlti, aice, i0,j0, puny)
      if len(iocn) == 0:
        # no ocn pnts:
        tsf_new = Tsfcn[:,j0,i0]*0.0 + Tfrz
      else:
        Tsf_adj = Tsfcn[:,jocn,iocn]
        tsf_new = np.nanmean(Tsf_adj, axis=1)
        tsf_new = np.where(tsf_new > Tsfc_max, Tsfc_max, tsf_new)
            
      vin_new = vin_old*0.0

    elif ai_old < puny and ai_new < puny:
      # Case: no ice --> no ice
      aice_case = "noice2noice"
      tsf_new = Tsfcn[:,j0,i0]
      vin_new = vin_old*0.0
    
    else:
      raise Exception(f"Unexpected case for ai_old={ai_old} and ai_new={ai_new}")  

    # Check that ice vol is correctly distributed across the ice thickn. cats:
    # If not - distribute across ice cats conserving aice and vice 
    # and matching ice cats
    ain_min = 1.e-8    # lower bound of ain(n) to avoid zeros
    vice_old = 0.      # ice vol / m2 _cell from old restart

    vice_old = np.sum(vin_old)
    vice_new = np.sum(vin_new)
    if ai_new > puny: 
      vice_m2ice_new = vice_new / ai_new  
    else:
      vice_m2ice_new = 0.

    vtot_target = vice_m2ice_fix * ai_new  # m3/m2_ice --> m3/m2_cell
    ain_new, vin_new = mc6util.adjust_thkncats_aice(ain_new, vin_new, vtot_target, \
                         hicat, dhi_min,  bnd_min=ain_min)
  
    ain_new = np.where(ain_new <= ain_min, 0., ain_new)
    vin_new = np.where(ain_new <= ain_min, 0., vin_new)

    if np.sum(ain_new) > (1.+1e-12): 
      logger.error(f"Ice conc > 1 {np.sum(ain_new)} j={j0} i={i0}")

    # Check ice cats:
    hin_new = np.divide(vin_new, ain_new, out=np.zeros_like(vin_new), where=ain_new != 0)
    cat_missed, hcat_new = mc6util.check_ithkn_cats(hicat, hin_new, ain_new)
    if cat_missed is not None:
      print(f"ipp={ipp} {aice_case} error in ice cats")
      print(f"vin_new={vin_new}, ain_new={ain_new}")
      raise Exception("Check ain_new, hin_new not in ice thkn cats")


    Tsfcn_new[:,j0,i0] = tsf_new
    vicen_new[:,j0,i0] = vin_new
    aicen_new[:,j0,i0] = ain_new

    # Update sice00?, qice00?
    for ilr in range(1, nilrs+1):
      varnum = f"{ilr:03d}"
      # Ice salinity by layers - compute S profile using BZ99 formulation:
      sice_lr = mc6util.sice_lr_cice4(ilr, nilrs, ain_new) 
      sice_old = sice[varnum][:,j0,i0]
      sice_new = np.where(sice_old < puny, sice_lr, sice_old)
      sice_new = np.where(ain_new < puny, 0., sice_new)

      sicen_new[varnum][:,j0,i0] = sice_new

      # ice enthalpy, use BL99
      # keep ice T below freezing T
      tice_max = sice_new*0 + Tfrz - 0.01
      qice_lr = mc6util.ice_enthalpy_BL99(tice_max, sice_new)
      qice_old = qice[varnum][:,j0,i0]
      qice_new = np.where(ain_old < puny, qice_lr, qice_old)
      qice_new = np.where(ain_new < puny, 0., qice_new)

      qicen_new[varnum][:,j0,i0] = qice_new

    # Check snow volume and ice_mean hsnow, add min snow if needed:
    vstot_old = np.sum(vsn_old)  # snow vol m3/m2_cell or cell mean snow thickn.
    if ai_new > puny:
      vsn_m2ice_new = vstot_old / ai_new  # mean snow thickn over ice
      if vsn_m2ice_new < hsnow_min:
        vstot_new = hsnow_min * ai_new
      elif vsn_m2ice_new > vsnow_m2ice_max:
        vstot_new = vsnow_m2ice_fix * ai_new
      else:
        vstot_new = vstot_old

      # Distribute evenly across cats:
      wts = ain_new/ai_new
      vsn_new = vstot_new * wts
    else:
      vsn_new = ain_new * 0.

    vsn_new = np.where(ain_new <= puny, 0., vsn_new)
    vsnon_new[:,j0,i0] = vsn_new

    # Update snow enthalpy: J/m3  
    # snow enthalpy should be: qsn_min <= qsn <= qsn_max
    # In theory, qsn_max = -rhos_Lfresh (latent heat of metling at 0C)
    # Make it a little lower to keep snow from melting right away
    #hsn_new = vsn_new / ain
    qsn = qsnon[:,j0,i0]        # enthalpy, J/kg < 0
    qsn_min = -rhos * Lfresh + (Tmin + 0.01) * cp_ice * rhos  # enth. of the coldest possible snow
    qsn_max = -rhos * Lfresh - 0.01 * cp_ice * rhos  # a little colder than 0C snow
    qsn_tsfc = -rhos * Lfresh + tsf_new * cp_ice * rhos # snow enth for surf. T
    qT0 = -Lfresh*rhos      # enth. of pure snow at 0C

    # Update new enthalpy of new snow:
    # Clip to min/max enthalpy, set to 0 where no snow:
    qsn_new = np.clip(qsn_tsfc, qsn_min, qsn_max)
    qsn_new = np.where(ain_new <= puny, 0., qsn_new)      # no ice
    qsn_new = np.where(vsn_new <= 0, 0., qsn_new)     # no snow
    qsnon_new[:,j0,i0] = qsn_new

  # Update data set:
  #ds_out = ds_out.assign(vsnon=vsnon_new, qsno001=qsnon_new)
  ds_out['vsnon'].values[:]   = vsnon_new
  ds_out['qsno001'].values[:] = qsnon_new
  ds_out['aicen'].values[:]   = aicen_new
  ds_out['vicen'].values[:]   = vicen_new
  ds_out['Tsfcn'].values[:]   = Tsfcn_new
  for ilr in range(1, nilrs+1):
    varnum = f"{ilr:03d}"
    ds_out[f"sice{varnum}"].values[:] = sicen_new[f"{varnum}"]
    ds_out[f"qice{varnum}"].values[:] = qicen_new[f"{varnum}"]

  # Sanity checking:
  assert ds_out["vsnon"].shape == vsnon_new.shape, "Check shape of vsnon "
  assert ds_out["qsno001"].shape == qsnon_new.shape, "Check shape of qsnon "

  check_print(aice, aicen_new, vicen_new, vsnon_new, ncat)

  return ds_out


def main():
  # Default values overridden by key arguments:
  vice_m2ice_max  = 10.             # max ice thickness (m3/m2_ice)
  vsnow_m2ice_max = 2.              # max snow depth (m3/m2_ice)  
  vice_m2ice_fix  = 1.              # new max ice thickness that replaces extreme ice thickness
  vsnow_m2ice_fix = 0.2             # new snow depth that replaces extreme snow thickness
  fyaml           = 'fix_cice.yaml' # yaml config file
  fsave           = 1               # =1: save new restart, =0: do not save new restart, debug mode

  parser = argparse.ArgumentParser()
  parser.add_argument("--use_yaml", action="store_true",
                      help="Use YAML values instead of script defaults for himax, hsmax, hifix, hsfix")
  parser.add_argument("--himax", type=float, 
                      help=f"Max ice thkn (m3/m2_ice), default={vice_m2ice_max:.2f}")
  parser.add_argument("--hsmax", type=float, 
                      help=f"Max snow depth (m3/m2_ice), default={vsnow_m2ice_max:.2f}")
  parser.add_argument("--hifix", type=float, 
                      help=f"New ice thkn (m3/m2_ice), default={vice_m2ice_fix:.2f}")
  parser.add_argument("--hsfix", type=float, 
                      help=f"New snow depth (m3/m2_ice), default={vsnow_m2ice_fix:.2f}")
  parser.add_argument("--fsave",type=int, help=f"Save new restart 0=no, 1=yes, default={fsave}",
                      choices=[0,1],default=fsave)
  parser.add_argument("--fyaml", help=f"Input yaml file with paths etc.", default=fyaml, type=str)
  parser.add_argument("--init",  help="Initialization date of the run YYYYMMDD", required=True, type=int)
  parser.add_argument("--ihr",   help="Initialization time, hour deault=0", default=0, type=int)
  parser.add_argument("--nmem",  help="Ensmb. member number", required=True, type=int)
  args = parser.parse_args()

  save_restart    = bool(args.fsave)
  fyaml           = args.fyaml
  init_hr         = args.ihr
  init_date       = args.init  if args.init is not None else None
  nmem            = args.nmem  if args.nmem is not None else None

  if not os.path.isfile(args.fyaml):
    raise FileNotFoundError(f"Missing YAML {args.fyaml}")

  cice_config = {}
  with open(args.fyaml) as ff:
    cice_config = safe_load(ff) or {}

  if args.use_yaml:
    vice_m2ice_max  = cice_config.get("max_hice", vice_m2ice_max)
    vsnow_m2ice_max = cice_config.get("max_hsnow", vsnow_m2ice_max)
    vice_m2ice_fix  = cice_config.get("fix_hice", vice_m2ice_fix)
    vsnow_m2ice_fix = cice_config.get("fix_hsnow", vsnow_m2ice_fix)
  else:
    vice_m2ice_max  = args.himax if args.himax is not None else vice_m2ice_max
    vsnow_m2ice_max = args.hsmax if args.hsmax is not None else vsnow_m2ice_max
    vice_m2ice_fix  = args.hifix if args.hifix is not None else vice_m2ice_fix
    vsnow_m2ice_fix = args.hsfix if args.hsfix is not None else vsnow_m2ice_fix

  pthrest     = cice_config["rest_path"].format(init_date=init_date, HH=init_hr, NR=nmem)
  pthrest_out = cice_config["out_path"].format(init_date=init_date, HH=init_hr, NR=nmem)

  # Get date of the 3hr back from the init_date:
  dnmb0 = mtime.rdate2datenum(init_date*100+init_hr)  # restart day nmb
  dnmbP = dnmb0 - 1./8. # assuming previous date/time is 3 hr back
  yrP, mmP, ddP, hrP = mtime.datevec(dnmbP, round_hrs=True)[:4]
 
  flrst_in  = f"{yrP}{mmP:02d}{ddP:02d}.{hrP:02d}0000.analysis.cice_model.res.nc"
  flrst_out = f"fix_{flrst_in}"  

  logger.info(
    f"Fixing CICE6 restart: hice: {vice_m2ice_max:.3f} --> {vice_m2ice_fix:.3f} "
    f"hsnow: {vsnow_m2ice_max:.3f} --> {vsnow_m2ice_fix:.3f}"
  )

  logger.info(f"restart input: {pthrest}/{flrst_in}")
  if save_restart:
    logger.info(f"restart output: {pthrest_out}/{flrst_out}")
  else:
    logger.info("Restart will not be saved")


  ds_out = None
  ds_out = fix_hsnow_ithkn(ds_out, vice_m2ice_max, vice_m2ice_fix, vsnow_m2ice_max, vsnow_m2ice_fix,
                           pthrest, flrst_in, logger)

  # Attributes:
  from datetime import datetime
  istep1_val = ds_out.attrs.get('istep1', None)
  ds_out.attrs.update({
      "title": f"CICE6 restart with fixed ice thkn and hsnow exceeding {vice_m2ice_max:.3f} m3/m2ice and {vsnow_m2ice_max:.3f} m3/m2ice",
      "source": "fix_cice6restart.py"
  })

  # Save:
  if save_restart:
    dflrst_out = os.path.join(pthrest_out, flrst_out)
    os.makedirs(pthrest_out, exist_ok=True)
    print(f"Saving CICE restart --> {dflrst_out}")
    ds_out.to_netcdf(dflrst_out, encoding={var: {'_FillValue': None} for var in ds_out.data_vars}, format='NETCDF3_64BIT')
    ds_out.close()

if __name__ == "__main__":
  main()


