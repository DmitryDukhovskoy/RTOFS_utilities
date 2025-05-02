"""
  Note SPEAR fields have now been interpoalted to NEP grid:
  calc_SPEAR_ice_clim.py


  Compute monthly SPEAR ice anomalies for ice conc and thickn. 
  inpoterlate onto NEP grid

  Saved months - f/cast months, not calendar
 
  usage: 
  interp_SPEAR_ice_anom_NEP.py --varnm iconc --YRS 2010 --YRE 2014 --MMI=1 --ensmb 1

  Only 12 months of the f/cast are saved in the SPEAR files
  mo - calendar month, depending on the init month, it may be at the end/start of the f/cast

"""
import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os
import importlib
import sys
import matplotlib.pyplot as plt
from yaml import safe_load
import argparse
import pickle

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
import mod_time as mtime
import mod_mom6 as mmom6
import mod_utils as mutil
import mod_colormaps as mclrmps
import mod_misc1 as mmisc
from mod_utils_fig import bottom_text
import mod_sis2_relax as msisrlx
importlib.reload(msisrlx)

# Default settings:
f_save = True
# Years in the relax file also used in the rlx file name:
YRS = 1993  # init yr
nyrs = 5  # numb. of years to derive climatologoies
YRE = YRS 
MINIT = [1,4,7,10]
ifld = 'iarea'  # ithkn, iarea
ENSMB = [x for x in range(1,11)]
varnm = 'iconc'
#ens_nmb = 1  # SPEAR ensemble #
plot_fields = False
f_save = True


parser = argparse.ArgumentParser()
parser.add_argument("--YRS", help="start year to calc anomalies: 1993, ..., 2020", type=int)
parser.add_argument("--YRE", help="end year to calc. anomalies: 1993, ..., 2020", type=int)
parser.add_argument("--MMI", help="init month of SPEAR f/cast: 1,4,7,10, default=all", type=int)
parser.add_argument("--varnm", help="field: ithkn or iarea/iconc", type=str)
args = parser.parse_args()

if args.YRS:
  YRS = args.YRS
if args.YRE:
  YRE = args.YRE
if args.MMI:
  MMI = args.MMI
  MINIT = [MMI]
if args.varnm:
  ifld = args.varnm
  if ifld=='iconc' or ifld=='iarea':
    ifld = 'siconc'
  elif ifld=='ithkn' or ifld=='ithk':
    ifld = 'sithick'

varnm = ifld

# Saved climatologies:
ICLIM=[[1993,1997],[1995,1999],[2000,2004],[2005,2009],[2010, 2014],[2015,2019],[2016,2020]]
ICLIM=np.array(ICLIM)

fconfig = 'config_nep.yaml'
with open(fconfig) as ff:
  config = safe_load(ff)
spear_dir = os.path.join(config['filesystem']['spear_month_ens'], 'monthly_clim')
# Check if mapping indices exist, gmapi:
dirgmapi = config['filesystem']['spear_mom_gmapi']
flgmaph  = 'spear2mom_NEP_full_gmapi_hpnt.nc'
dflgmaph = os.path.join(dirgmapi, flgmaph)
# h-point indices
if not os.path.isfile(dflgmaph):
  print(f'Mapping indices hpnt are missing, {dflgmaph}, run find_SPEAR_NEP_gmapi.py ...')
  raise Exception('Quitting ...')
dsh = xarray.open_dataset(dflgmaph)
IMOM = dsh['indx_nep'].data
JMOM = dsh['jndx_nep'].data
INDX = dsh['indx_spear'].data
JNDX = dsh['jndx_spear'].data
LON  = dsh['lonh_spear'].data
LAT  = dsh['lath_spear'].data

# NEP grid:
fyaml = 'paths_seasfcst.yaml'
with open(fyaml) as ff:
  pthseas = safe_load(ff)

expt_name = "seasonal_daily"
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
jdim, idim = HH.shape

LMsk = np.where(HH>=0, 0, 1)
# Mask out southern lats:
LMsk[:551,:] = 0
LMsk[:,:47] = 0


# Read SPEAR relax. fields:
icc = 0
for YRI in range(YRS,YRE+1):
  for MMI in MINIT:
    AENSMB = None 
    for ens_nmb in ENSMB:  
      print(f'Processing {YRI}')
      pthdata = f'/work/Dmitry.Dukhovskoy/tmp/spear_subset/{YRI}/ens{ens_nmb:02d}'
      flthkn = f'NEP_spear_{YRI}{MMI:02d}.sithick.nc'
      dfthkn = os.path.join(pthdata,flthkn)
      dspear_thkn = xarray.open_dataset(dfthkn)
      flconc = f'NEP_spear_{YRI}{MMI:02d}.siconc.nc'
      dfconc = os.path.join(pthdata,flconc)
      dspear_conc = xarray.open_dataset(dfconc)

      # Calendar months of the f/cast period:
      mcal_fcast = np.arange(MMI,MMI+12)
      mcal_fcast = np.where(mcal_fcast>12, mcal_fcast-12, mcal_fcast)

      # Get climatology not interpolated to NEP:
      iclm = np.where((ICLIM[:,0] <= YRI) & (ICLIM[:,1] >= YRI))[0]
      assert(len(iclm)>0), f'Could not find clim. time window for {YRI}'
      iclm = iclm[0]
      YRC1 = ICLIM[iclm,0]
      YRC2 = ICLIM[iclm,1] 
      pthpkl = '/work/Dmitry.Dukhovskoy/anls_output/spear_ice'
      flclim = f'spear_siconc_clim_{YRC1}_{YRC2}_MI{MMI:02d}e{ens_nmb:02d}.nc'
      dflclim = os.path.join(pthpkl,flclim)
      print(f'Getting climtology --> {dflclim}')
      # Note SPEAR clim. starts from month=MMI !!!
      ds_clim = xarray.open_dataset(dflclim)
      Aclim = ds_clim['siconc'].data
      LON   = ds_clim['longitudes'].data
      LAT   = ds_clim['latitudes'].data
      mcal_clim  = ds_clim['calend_months'].data  # clim months, should match clim_f/cast

      for itime in range(12):
        MMF = itime+1  # forecast lead time in months, 1,..., 12
        mm_cal = mcal_fcast[itime]
        C2d = dspear_conc['siconc'].isel(time=itime).data
        iclim = np.where(mcal_clim == mm_cal)[0][0]
        Anom = C2d - Aclim[iclim,:,:].squeeze()
        # Put anomalies in calendar months order:
        ical = mm_cal-1
        if AENSMB is None:
          dim1, dim2 = C2d.shape
          AENSMB = np.zeros((12,dim1,dim2))
        AENSMB[ical,:,:] = AENSMB[ical,:,:] + Anom
  
    nensmb = len(ENSMB)
    AENSMB = AENSMB/nensmb

    # Interpolate ensemble mean anomalies
    # Arrange in calendar months order
    Ianom = np.zeros((12,jdim,idim))
    for ical in range(12):
      mm_cal = ical+1
      mm_lead = msisrlx.cal_mo_to_fcast(MMI,mm_cal)
      print(f'Interpolating anomalies, {YRI} MMI={MMI} cal_month={mm_cal}, ' +\
            f'fcast_lead_month={mm_lead}')
      Anom = AENSMB[ical,:,:]
      amin1 = np.nanmin(Anom)
      amax1 = np.nanmax(Anom)

      # Interpolate to NEP grid:
      Anomi = msisrlx.interp2Dfld(Anom, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)
      Anomi = np.where(np.isnan(Anomi), 0., Anomi)
      Anomi = np.where(HH>=0, np.nan, Anomi)
      Ianom[ical,:,:] = Anomi
      amin2 = np.nanmin(Anomi)
      amax2 = np.nanmax(Anomi)

      print(f'min/max before interp: {amin1:.3f}/{amax1:.3f}, after: {amin2:.3f}/{amax2:.3f}') 

    mcal = np.arange(1,13)
    kdim,jdim,idim = Ianom.shape
    darr_ice = xarray.DataArray(Ianom, dims=("months","jdim","idim"), \
                      coords={"months": mcal, \
                              "jdim": np.arange(jdim), \
                              "idim": np.arange(idim)})
    darr_months = xarray.DataArray(mcal, dims=("months"), \
                         coords={"months": mcal})
    dset = xarray.Dataset({"calend_months": darr_months, f"{varnm}_anom": darr_ice})

    if f_save:
      flanom = f'spear_{varnm}_monanom_ensmean_{YRI}{MMI:02d}.nc'
      dflanom = os.path.join(pthpkl,flanom)
      print(f'Dumping {varnm} anomalies --> {dflanom}')
      dset.to_netcdf(dflanom, format='NETCDF3_64BIT', engine='netcdf4')


f_check = False
if f_check:
  CLRS = [[0.6, 0.02, 0.6],
          [0.2, 0.38, 1],
          [0., 0.8, 0.5],
          [0.2,1.,0.8],
          [1, 1, 1],
          [1, 0.9, 0.85],
          [1, 0.4, 0.4],
          [0.9, 0.6,0],
          [0.6, 0.2, 0]]

  tcmp = mclrmps.colormap_posneg_uneven(CLRS)
  tmin = -2.
  tmax = 5.

  clrmp_dlt = mclrmps.colormap_ssh(cpos='YlOrRd', cneg='PuBuGn_r')
  clrmp_dlt.set_bad(color=[0.6,0.6,0.6])
  rmin = -1.
  rmax = 1.


  clrmp.set_bad(color=[0.2, 0.2, 0.2])

  plt.ion()

  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  img = ax1.pcolormesh(Anomi, cmap=clrmp_dlt, vmin=rmin, vmax=rmax)

  sttl = f'SPEAR {varnm} anom wrt {YRC1}-{YRC2} mean'
  ax1.set_title(sttl)

  ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                     0.02, ax1.get_position().height])
  # extend: min, max, both
  clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
  ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
  #  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.1f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=12)

  btx = 'interp_SPEAR_ice_anom_NEP.py'
  bottom_text(btx, pos=[0.2, 0.01])


