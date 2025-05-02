"""
  Compute monthly SPEAR ice climatologies for ice conc and thickn. 
  for all esnembles and init. dates to derive ensemble-mean ice anomalies etc.  

  Extract and subsample for NEP domain using bash script:
  ./subset_spear_ice.sh 2010 2010 1 1

  Note that year 2010 and month 1 designate the initialization time
  SPEAR monthly fields have 12 f/cast months in each file

  Do not interpolate to NEP, to save time
  compute anomalies on SPEAR grid, derive ensmble mean then interpolate
  Default: interp =0

  usage: 
  calc_SPEARensmb_ice_clim.py --varnm iconc --YRS 2010 --nyears 5 --MMI=1 --ensmb 1 --interp 1

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
nyrs = 5
YRE = YRS + (nyrs-1)
MINIT = [1,4,7,10]
#MMI = 1     # init month
ifld = 'iarea'  # ithkn, iarea
ENSMB = [x for x in range(1,11)]
#ens_nmb = 1  # SPEAR ensemble #
interp_NEP = False # interpolate to NEP and save on both grids (SPEAR & NEP), otherwise - only SPEAR

parser = argparse.ArgumentParser()
parser.add_argument("--YRS", help="start year for deriving climat.: 1993, ..., 2022", type=int)
parser.add_argument("--nyears", help=f"number of years to derive climat., default={nyrs}", type=int)
parser.add_argument("--MMI", help="init month of SPEAR f/cast: 1,4,7,10, default=all", type=int)
parser.add_argument("--ensmb", help="ensemble number, 1,..., 10, default=all", type=int)
parser.add_argument("--varnm", help="field: ithkn or iarea/iconc", type=str)
parser.add_argument("--interp", help="=1: interp to NEP grid and save on both NEP and SPEAR, default=0", type=int) 
args = parser.parse_args()

if args.YRS:
  YRS = args.YRS
  YRE = YRS + (nyrs-1)
if args.nyears:
  nyrs = args.nyears
  YRE = YRS + (nyrs-1)
if args.MMI:
  MMI = args.MMI
  MINIT = [MMI]
if args.varnm:
  ifld = args.varnm
  if ifld=='iconc' or ifld=='iarea':
    ifld = 'siconc'
  elif ifld=='ithkn' or ifld=='ithk':
    ifld = 'sithick'
if args.ensmb:
  ens_nmb = args.ensmb
  ENSMB = [ens_nmb]
if args.interp:
  if args.interp == 1:
    interp_NEP = True
  else:
    interp_NEP = False

varnm = ifld

if not f_save:
  print(f'WARNING: fields wont be saved, f_save flag is off !!!\n')

# Read SPEAR fields:
# Note that SPEAR thickness is not "volume/unit area" (like in CICE output or PIOMAS)
# and needs to be multiplied by partial area
for ens_nmb in ENSMB:
  for MMI in MINIT: 
    icc = 0
    for YRI in range(YRS,YRE+1):
      print(f'Processing {YRI} MMI{MMI:02d}e{ens_nmb:02d}')
      pthdata = f'/work/Dmitry.Dukhovskoy/tmp/spear_subset/{YRI}/ens{ens_nmb:02d}'
      flthkn = f'NEP_spear_{YRI}{MMI:02d}.sithick.nc'
      dfthkn = os.path.join(pthdata,flthkn)
      dspear_thkn = xarray.open_dataset(dfthkn)
      flconc = f'NEP_spear_{YRI}{MMI:02d}.siconc.nc'
      dfconc = os.path.join(pthdata,flconc)
      dspear_conc = xarray.open_dataset(dfconc)
      # Assumed that rec #1 = init month
      #Time = ds_spear['time'].data
      #TM = mmisc.convert_nptime_to_datenum(Time)
      #dnmb0 = mtime.datenum([YR0,MM0,15,12])
      mcal = np.arange(MMI,MMI+12)
      mcal = np.where(mcal>12, mcal-12, mcal)
      #D = abs(mcal-MM0)
      #itime = np.argmin(D)

      if icc == 0:
        xh = dspear_thkn['xh'].data
        yh = dspear_thkn['yh'].data
        xdim = len(xh)
        ydim = len(yh)
        ASUM = np.zeros((12,ydim,xdim))
        
      # Process by f/casts months
      for itime in range(12):
        H2d = dspear_thkn['sithick'].isel(time=itime).data
        C2d = dspear_conc['siconc'].isel(time=itime).data
        if ifld == 'siconc':
          A2d = C2d.copy()
        elif ifld == 'sithick':
          A2d = H2d*C2d         # ice m ---> m3/m2 

        hmax = np.nanmax(A2d)
        hmin = np.nanmin(A2d)
        print(f'N={itime+1} M={mcal[itime]} {ifld} min/max: {hmin:.2f}/{hmax:.2f}')

        ASUM[itime,:,:] = ASUM[itime,:,:] + A2d
        if icc == 0:
          LON = dspear_conc['GEOLON'].data
          LAT = dspear_conc['GEOLAT'].data

      icc += 1

    ASUM = ASUM / icc

    if f_save:
      dim1,dim2,dim3 = ASUM.shape
      darr_ice = xarray.DataArray(ASUM, dims=("months","jdim","idim"), \
                        coords={"months": np.arange(dim1), \
                                "jdim": np.arange(dim2), \
                                "idim": np.arange(dim3)})
      darr_lon = xarray.DataArray(LON, dims=("jdim","idim"),\
                        coords={"jdim": np.arange(dim2), \
                                "idim": np.arange(dim3)})
      darr_lat = xarray.DataArray(LAT, dims=("jdim","idim"),\
                        coords={"jdim": np.arange(dim2), \
                                "idim": np.arange(dim3)})
      darr_months = xarray.DataArray(mcal, dims=("months"), \
                           coords={"months": np.arange(dim1)})
      dset = xarray.Dataset({"calend_months": darr_months, \
                             "longitudes": darr_lon, \
                             "latitudes": darr_lat, \
                             f"{varnm}": darr_ice})
      dset['calend_months'].attrs['long_name']='Calendar months during the forecast'
      if ifld == 'siconc':
        dset[f'{varnm}'].attrs['long_name']='ice partial area'
      elif ifld == 'sithick':
        dset[f'{varnm}'].attrs['long_name']='mean thickness or ice volume per unit area, m3/m2'
        
      # Add global attributes:
      dset.attrs.update({
        "info": f"SPEAR {varnm} climatology on SPEAR grid subset of the N.Pacfic region",
        "code": "calc_SPEARensmb_ice_clim.py"
      })
        
      pthpkl = '/work/Dmitry.Dukhovskoy/anls_output/spear_ice'
      flint  = f'spear_{varnm}_clim_{YRS}_{YRE}_MI{MMI:02d}e{ens_nmb:02d}.nc'
      dflint = os.path.join(pthpkl,flint)
      print(f'Dumping SPEAR not-interpolated {varnm} climtology --> {dflint}')
      dset.to_netcdf(dflint, format='NETCDF3_64BIT', engine='netcdf4')

    if interp_NEP:
      print("Interpolating to NEP grid ...")

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

      AINT = np.zeros((12,jdim,idim))
      fconfig = 'config_nep.yaml'
      with open(fconfig) as ff:
        config = safe_load(ff)
      # Check if mapping indices exist, gmapi, should be created before interpolation:
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
      LON2  = dsh['lonh_spear'].data
      LAT2  = dsh['lath_spear'].data

      # Make sure that gmapi are for this SPEAR subset:
      D = np.max(np.abs(LON-LON2) + np.abs(LAT-LAT2))
      if D > 1.e-6:
        raise Exception ('Saved gmapi may not be for this SPEAR LON/LAT, check subset regions ...')
       
      for imo in range(12):
        print(f"  month {imo}")
        # Interpolate to NEP grid:
        A2d = ASUM[imo,:,:].squeeze()
        A2di = msisrlx.interp2Dfld(A2d, IMOM, JMOM, INDX, JNDX, LMsk, LON, LAT, hlon, hlat)
        A2di = np.where(np.isnan(A2di), 0., A2di)
        A2di = np.where(HH>=0, np.nan, A2di)

        AINT[imo,:,:] = A2di

      kdim,jdim,idim = AINT.shape
      darr_ice = xarray.DataArray(AINT, dims=("months","jdim","idim"), \
                        coords={"months": np.arange(kdim), \
                                "jdim": np.arange(jdim), \
                                "idim": np.arange(idim)})
      darr_months = xarray.DataArray(mcal, dims=("months"), \
                           coords={"months": np.arange(kdim)})
      dset = xarray.Dataset({"calend_months": darr_months, f"{varnm}": darr_ice})
      dset['calend_months'].attrs['long_name']='Calendar months during the forecast'
      if ifld == 'siconc':
        dset[f'{varnm}'].attrs['long_name']='ice partial area'
      elif ifld == 'sithick':
        dset[f'{varnm}'].attrs['long_name']='mean cell thickness or ice volume per unit area, m3/m2'

      # Add global attributes:
      dset.attrs.update({
        "info": f"Interpolated SPEAR {varnm} climatology to NEP SIS2 grid",
        "code": "calc_SPEAR_ice_clim.py"
      })

      if f_save:
        flint  = f'spear_interpNEP_{varnm}_clim_{YRS}_{YRE}_MI{MMI:02d}e{ens_nmb:02d}.nc'
        dflint = os.path.join(pthpkl,flint)
        print(f'Dumping interpolated {varnm} climtology --> {dflint}')
        dset.to_netcdf(dflint, format='NETCDF3_64BIT', engine='netcdf4')
     

match ifld:
  case('sithick'):
    clrmp = mclrmps.colormap_ice_thkn()
    rmin = 0.
    rmax = 4.
  case('siconc'):
    clrmp = mclrmps.colormap_conc()
    rmin = 0.
    rmax = 1.

clrmp.set_bad(color=[0.2, 0.2, 0.2])


def plot_ice(fgnmb, m, xR, yR, A2d, clrmp, rmin, rmax, sttl):
  fig1 = plt.figure(fgnmb,figsize=(9,8))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  m.drawcoastlines()
  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

  img = ax1.pcolormesh(xR, yR, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  #  img = ax1.pcolormesh(RLXHR, cmap=clrmp)

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

  btx = 'calc_SPEAR_ice_clim.py'
  bottom_text(btx, pos=[0.2, 0.01])


# Check
f_check = False

if f_check:
  plt.ion()

  MM0 = 12  # month to check
  D = abs(mcal-MM0)
  it0 = np.argmin(D)
  A2d = ASUM[it0,:,:].squeeze()
  #A2d = np.where(HH>=0, np.nan, A2d)

  # Stereographic Map projection:
  from mpl_toolkits.basemap import Basemap, cm
  #m = Basemap(width=5000*1.e3,height=5000*1.e3, resolution='l',\
  #            projection='stere', lat_ts=50, lat_0=62, lon_0=-165)
  m = Basemap(width=3300*1.e3,height=3700*1.e3, resolution='l',\
            projection='stere', lat_ts=60, lat_0=65, lon_0=-175)

  xR, yR = m(LON, LAT)

  fgnmb=1
  sttlS = f'SPEAR clim init: {YRI}/{MMI:02d}-e{ens_nmb:02d}, {ifld} {YRS}-{YRE} {MM0}'
  plot_ice(fgnmb, m, xR, yR, A2d, clrmp, rmin, rmax, sttlS)

