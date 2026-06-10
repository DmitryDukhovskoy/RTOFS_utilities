"""
  Check snow depth at a grid point

  From a debug run with saved variables in ice thickness categories

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
import matplotlib.colors as colors
from yaml import safe_load
from mpl_toolkits.basemap import Basemap, cm
import argparse

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

from mod_utils_fig import bottom_text
import mod_time as mtime
import mod_colormaps as mclrmps
import mod_mom6 as mmom6
import mod_misc1 as mmisc

init_date = 20250701
init_hr = 0

# 61 - debug run
enmb = 61

parser = argparse.ArgumentParser()
parser.add_argument("--init", help=f"init date", choices=[20240701, 20250701], default=init_date, type=int)
parser.add_argument("--ihr", help=f"init hour, default {init_hr}", type=int)
parser.add_argument("--fday1", help=f"forecast day to plot: 0, 1, 2, ... , =0 - init. cond.", type=int, required=True)
parser.add_argument("--fday2", help=f"forecast day to plot: 0, 1, 2, ... , =0 - init. cond.", type=int, required=True)
parser.add_argument("--enmb", help="experiment number: 0, 1, 2, ...", default=enmb, type=int)
args = parser.parse_args()

init_date = args.init if args.init else init_date
init_hr   = args.ihr if args.ihr else init_hr
fday1     = args.fday1
fday2     = args.fday2
enmb      = args.enmb
rho_snow = 300.

def path_names(fday, init_date, init_hr, enmb):
  plot_init = fday == 0  # initial conditions
  
  dnmbI = mtime.rdate2datenum(init_date*100+init_hr)  # init. day nmb
  YRI,MMI,DDI,hrI = mtime.datevec(dnmbI, round_hrs=True)[:4]

  if plot_init:
    dnmb0 = dnmbI
  else:
    dnmb0 = dnmbI + fday-1                           # day to plot

  yr0,mm0,dd0,hr0 = mtime.datevec(dnmb0, round_hrs=True)[:4]
  nsec0 = int(hr0 * 3600)
  YR,MM,DD = mtime.datevec(dnmb0)[:3]

  pthoutp = f"/gpfs/f6/sfs-cpu/scratch/Dmitry.Dukhovskoy/sfs_C192mx025_cice_test/expt{enmb:02d}/cice6"
  if plot_init:
    flinp = f"iceh_ic.{yr0}-{mm0:02d}-{dd0:02d}-{nsec0:05d}.nc"
  else:
    flinp = f"iceh.{yr0}-{mm0:02d}-{dd0:02d}.nc"

  dflice = os.path.join(pthoutp,flinp)
  #print(f"Reading {YR}/{MM}/{DD}, SFS init {YRI}/{MMI:02d}/{DDI:02d}\n  {dflice}")

  return dflice

def read_cice(fday, HH, init_date, init_hr, enmb):
  dflice = path_names(fday, init_date, init_hr, enmb)
  print(f"Reading {dflice}")

  with xarray.open_dataset(dflice) as dcice:
    vsno_m2c = dcice['hs_d'].data.squeeze()   # grid cell mean snow vol, m3/m2_cell
    aicen = dcice['aicen_d'].values.squeeze()  # ice area in cats
    vicen = dcice['vicen_d'].values.squeeze()  # ice vol m2_ice in cats
    vsnon = dcice['vsnon_d'].values.squeeze()  # snow vol m2_snow (!) in cats
    snfrn = dcice['snowfracn_d'].values.squeeze() # snow frac m2_cell, in cats

  vsno_m2c[HH >= 0] = np.nan

  return vsno_m2c, aicen, vicen, vsnon, snfrn

def read_cicefld_gridpnt(fday, init_date, init_hr, varnm, jj0, ii0, enmb):
  dflice = path_names(fday, init_date, init_hr, enmb)
  #print(f"Reading {varnm} from {dflice}")

  with xarray.open_dataset(dflice) as dcice:
    FLD = dcice[varnm].values.squeeze()

  if FLD.ndim == 2:
    apnt = FLD[jj0,ii0]
  elif FLD.ndim == 3:
    apnt = FLD[:,jj0,ii0]
  elif FLD.ndim == 4:
    apnt = FLD[:,:,jj0,ii0]
  else:
    raise ValueError(
        f"Unsupported dimensions for {varnm}: {FLD.shape}"
    )

  return apnt

def read_cicefld_attr(fday, init_date, init_hr, varnm, jj0, ii0, enmb):
  # This is for scalar output only, i.e. not for (ncat,:,:) fields
  dflice = path_names(fday, init_date, init_hr, enmb)
  #print(f"Reading {varnm} from {dflice}")

  with xarray.open_dataset(dflice) as dcice:
    if varnm not in dcice:
      raise KeyError(f"{varnm} not found in {dflice}")
    FLD = dcice[varnm].values.squeeze()
    units = dcice[varnm].attrs.get('units', '')
    long_name = dcice[varnm].attrs.get('long_name', '')

  if FLD.ndim == 2:
    apnt = FLD[jj0,ii0]
  elif FLD.ndim == 3:
    apnt = FLD[:,jj0,ii0]
  else:
    raise ValueError(
        f"Unsupported dimensions for {varnm}: {FLD.shape}"
    )

  return apnt, units, long_name

def print_diagn_pnt(fday1, fday2, init_date, init_hr, jj0, ii0, enmb, VARS):
  lname_w = 40
  units_w = 10
  varnm_w = 15

  def fmt_val(x):
    if 0 < abs(x) < 1.e-6:
      return f"{x:12.5e}"
    else:
      return f"{x:12.7f}"

  # Print 1 pnt diagn for 2D flds (not cats)
  sfx = '_d' # daily average fields
  for var in VARS:
    varnm = f"{var}{sfx}"

    val1, units, lname = read_cicefld_attr(fday1, init_date, init_hr, varnm, jj0, ii0, enmb)
    val2, _, _ = read_cicefld_attr(fday2, init_date, init_hr, varnm, jj0, ii0, enmb)
    
    #print(f"Reading {varnm}: {units}  {lname}")
    # Truncate if long name is too long: 
    print(
        f"{lname[:lname_w]:<{lname_w}} "
        f"{units:<{units_w}} "
        f"{var:<{varnm_w}} :  "
        f"{fmt_val(val1)}   {fmt_val(val2)}"
    )

  return
  
syst_info = os.uname() 
machine = syst_info.nodename
  
if 'dtn' in machine:
  print("Running on DTN node:", machine)
  node_nm = "dtn"
elif 'gaea' in machine:
  print("Running on Gaea compute node:", machine)
  node_nm = "gaea"
else:
  print("Unknown machine:", machine)

fyaml = 'paths_ufs.yaml'
with open(fyaml) as ff:
  pths_ufs = safe_load(ff)

# Get MOM6 grid
pthgrid = pths_ufs[node_nm]["MOM6"]["pthgrid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

hlon, hlat = mmom6.read_mom6grid(dfgrid_mom, grdpnt='hgrid')


with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape

Vsno1_m2c, Aicen1, Vicen1, Vsnon1, Snfrn1 = read_cice(fday1, HH, init_date, init_hr, enmb)
Vsno2_m2c, Aicen2, Vicen2, Vsnon2, Snfrn2 = read_cice(fday2, HH, init_date, init_hr, enmb)

clrmp = mclrmps.colormap_temp()
rmin = 0.
rmax = 0.1

clrmp.set_bad(color=[0.1, 0.1, 0.1])
cntr_clr = [0.9,0.,1]

plt.ion()

print("Plotting ...")


# snowfracn can be > 1 ???
# np.nanmax(Snfrn1) = 2. ?
# Plot aggregated vsno(n)*snowfracn(n) - should match Vsno?_m2c 
# but it does not 
# To match, need to divide Snfrn1 / 2, then the 2 fields match - IC field
# Somehow this does not work for not IC field. 
aggr_vsn1 = np.sum(Vsnon1 * Snfrn1/2, axis=0).squeeze()
aggr_vsn2 = np.sum(Vsnon2 * Snfrn2, axis=0).squeeze()

A2d = Vsno2_m2c  # grid cell mean snow depth (vol snow / m2_grid)
#A2d = aggr_vsn2



sttl = f"SFS expt{enmb} init {init_date} FDAY={fday2}, grid cell mean hsnow"

fig1 = plt.figure(1,figsize=(9,9))
plt.clf()
ax1 = plt.axes([0.08, 0.1, 0.8, 0.8])

img = ax1.pcolormesh(A2d, cmap=clrmp, vmin=rmin, vmax=rmax, shading='auto')

ax1.set_aspect('equal', adjustable='box')
ax1.set_ylim(800, 1080)
ax1.set_xlim(100, 670)
ax1.set_title(sttl)

# Grid point analysis
# Rapid snow change:
ii0 = 303
jj0 = 1062

# Low snow change:
#ii0 = 311
#jj0 = 1070

ax1.plot(ii0, jj0, marker='o', color=[1, 0.0, 0.8])

ax3 = fig1.add_axes([0.2, 0.1, 0.6, 0.02])
clb = plt.colorbar(img, cax=ax3, orientation='horizontal', extend='max')
ax3.xaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax3.set_xticklabels(ax3.get_xticks())
ticklabs = clb.ax.get_xticklabels()
clb.ax.set_xticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=12)
clb.ax.tick_params(direction='in', length=12)


btx = 'debug_hsnow.py'
bottom_text(btx)

# Convert rate of snow loss (liquid) to ~ cm of snow (300 kg/m3) over 1 day:
# Rate of snow loss to leads (kg / (m2*sec))
fsloss1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "fsloss_d", jj0, ii0, enmb)
fsloss2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "fsloss_d", jj0, ii0, enmb)

cff = 3600*24*100/rho_snow  # cm of snow / day
fsloss1_cm = fsloss1 * cff
fsloss2_cm = fsloss2 * cff


# Get values in the grid
hs1 = Vsno1_m2c[jj0,ii0]  # Grid cell mean snow thickness = m3(snow) / m2_cell
hs2 = Vsno2_m2c[jj0,ii0]
aicen1 = Aicen1[:,jj0,ii0]
aicen2 = Aicen2[:,jj0,ii0]
vicen1 = Vicen1[:,jj0,ii0]
vicen2 = Vicen2[:,jj0,ii0]
vsnon1 = Vsnon1[:,jj0,ii0]
vsnon2 = Vsnon2[:,jj0,ii0]
snfrn1 = Snfrn1[:,jj0,ii0]  # <-- not sure, why it is > 1 ???
snfrn2 = Snfrn2[:,jj0,ii0]  # <-- ???

# snow int. T:
tsnz1 = read_cicefld_gridpnt(fday1, init_date, init_hr, 'Tsnz_d', jj0, ii0, enmb) 
tsnz2 = read_cicefld_gridpnt(fday2, init_date, init_hr, 'Tsnz_d', jj0, ii0, enmb) 
# Liquid mass in snow:
snliq1 = read_cicefld_gridpnt(fday1, init_date, init_hr, 'smassliqn_d', jj0, ii0, enmb) 
snliq2 = read_cicefld_gridpnt(fday2, init_date, init_hr, 'smassliqn_d', jj0, ii0, enmb) 

# 4D fields
# Ice int. T:
tinz1 = read_cicefld_gridpnt(fday1, init_date, init_hr, 'Tinz_d', jj0, ii0, enmb)
tinz2 = read_cicefld_gridpnt(fday2, init_date, init_hr, 'Tinz_d', jj0, ii0, enmb)

aice1 = np.sum(aicen1)
aice2 = np.sum(aicen2)
hicen1 = vicen1 / aicen1
hicen2 = vicen2 / aicen2
vsno1 = np.sum(vsnon1 * aicen1)  # m3/m2_ice
vsno2 = np.sum(vsnon2 * aicen2)

VARS = ['fsloss', 'hi', 'hs', 'snowfrac', 'Tair', 'Tsfc', 'sice', 
        'fswdn', 'fswup','flwdn', 'flwup', 'snow', 'rain', 'frzmlt', 'albsni',
        'albsno', 'albice', 'albpnd', 'flat', 'fsens', 'fsurf_ai', 'congel',
        'frazil', 'snoice', 'dsnow', 'melts', 'meltt',
        'meltb', 'meltl', 'fbot', 'fhocn', 'fswthru', 'apond',
        'hpond', 'ipond',  'fsloss','rsnw']

print(f" ------------------------")
print(f"Grid pnt: j={jj0} i={ii0}\n HH={HH[jj0,ii0]:.2f}, lon={hlon[jj0,ii0]:.2f} lat={hlat[jj0,ii0]:.2f}")
print_diagn_pnt(fday1, fday2, init_date, init_hr, jj0, ii0, enmb, VARS)

# 3D fields with ice cats.
str = f"vsnon j={jj0} i={ii0}:"
mmisc.print_cols(vsnon1, vsnon2, str=str)

str = f"aicen j={jj0} i={ii0}:"
mmisc.print_cols(aicen1, aicen2, str=str)

str = f"vicen j={jj0} i={ii0}:"
mmisc.print_cols(vicen1, vicen2, str=str)

str = f"hice j={jj0} i={ii0}:"
mmisc.print_cols(hicen1, hicen2, str=str)

str = f"tsnz snow int T j={jj0} i={ii0}:"
mmisc.print_cols(tsnz1, tsnz2, str=str)

# 4D fields (cats, layers, j ,i)
ncat, nilrs = tinz1.shape
print(f"tinz, Ice int T by cats:")
for icat in range(1,ncat+1):
  str = f'cat = {icat}'
  mmisc.print_cols(tinz1[icat-1,:], tinz2[icat-1,:], str=str)


#   =======
# Get fields:
# Rate of snow loss to lead, kg / (m2*sec)
fsloss1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "fsloss_d", jj0, ii0, enmb)
fsloss2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "fsloss_d", jj0, ii0, enmb)
# Grid cell mean ice thickness = m3(ice) / m2_cell
hi1     = read_cicefld_gridpnt(fday1, init_date, init_hr, "hi_d", jj0, ii0, enmb)
hi2     = read_cicefld_gridpnt(fday2, init_date, init_hr, "hi_d", jj0, ii0, enmb)
#Grid cell mean snow fraction:
snfr1  = read_cicefld_gridpnt(fday1, init_date, init_hr, "snowfrac_d", jj0, ii0, enmb)
snfr2  = read_cicefld_gridpnt(fday2, init_date, init_hr, "snowfrac_d", jj0, ii0, enmb)
# Sn/ice Tsfc:
tsfc1  = read_cicefld_gridpnt(fday1, init_date, init_hr, "Tsfc_d", jj0, ii0, enmb)
tsfc2  = read_cicefld_gridpnt(fday2, init_date, init_hr, "Tsfc_d", jj0, ii0, enmb)
# Bulk ice S:
sice1  = read_cicefld_gridpnt(fday1, init_date, init_hr, "sice_d", jj0, ii0, enmb)
sice2  = read_cicefld_gridpnt(fday2, init_date, init_hr, "sice_d", jj0, ii0, enmb)
# Down sh/wave:
fswdn1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "fswdn_d", jj0, ii0, enmb)
fswdn2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "fswdn_d", jj0, ii0, enmb)
# Down l/wave:
flwdn1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "flwdn_d", jj0, ii0, enmb)
flwdn2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "flwdn_d", jj0, ii0, enmb)
# Snofall (liq. water. equiv) --> snow:
snfall1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "snow_d", jj0, ii0, enmb) * 1000/rho_snow
snfall2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "snow_d", jj0, ii0, enmb) * 1000/rho_snow
# Frz-melt pot. 
frzmlt1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "frzmlt_d", jj0, ii0, enmb)
frzmlt2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "frzmlt_d", jj0, ii0, enmb)
# Snow/ice broad-band albedo:
albsni1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "albsni_d", jj0, ii0, enmb)
albsni2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "albsni_d", jj0, ii0, enmb)
# Snow albedo:
albsno1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "albsno_d", jj0, ii0, enmb)
albsno2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "albsno_d", jj0, ii0, enmb)
# Bare ice albedo:
albice1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "albice_d", jj0, ii0, enmb)
albice2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "albice_d", jj0, ii0, enmb)
# Laten heat flux:
flat1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "flat_d", jj0, ii0, enmb)
flat2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "flat_d", jj0, ii0, enmb)
# Sensible heat flux:
fsens1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "fsens_d", jj0, ii0, enmb)
fsens2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "fsens_d", jj0, ii0, enmb)
# Congel ice growth:
congel1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "congel_d", jj0, ii0, enmb)
congel2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "congel_d", jj0, ii0, enmb)
# Frazil ice growth:
frazil1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "frazil_d", jj0, ii0, enmb)
frazil2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "frazil_d", jj0, ii0, enmb)
# Snow-ice formation (flooded snow --> ice):
snoice1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "snoice_d", jj0, ii0, enmb)
snoice2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "snoice_d", jj0, ii0, enmb)
# Snow formation
dsnow1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "dsnow_d", jj0, ii0, enmb)
dsnow2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "dsnow_d", jj0, ii0, enmb)
# Top snow melt:
melts1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "melts_d", jj0, ii0, enmb)
melts2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "melts_d", jj0, ii0, enmb)
# Top ice melt
meltt1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "meltt_d", jj0, ii0, enmb)
meltt2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "meltt_d", jj0, ii0, enmb)
# Basal ice melt:
meltb1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "meltb_d", jj0, ii0, enmb)
meltb2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "meltb_d", jj0, ii0, enmb)
# Lateral ice melt
meltl1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "meltl_d", jj0, ii0, enmb)
meltl2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "meltl_d", jj0, ii0, enmb)
# Heat flux ice to ocean (cpl):
fhocn1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "fhocn_d", jj0, ii0, enmb)
fhocn2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "fhocn_d", jj0, ii0, enmb)
# Sh/wave through ice->ocean:
fswthru1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "fswthru_d", jj0, ii0, enmb)
fswthru2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "fswthru_d", jj0, ii0, enmb)
# Melt pond fraction of sea ice:
apond1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "apond_d", jj0, ii0, enmb)
apond2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "apond_d", jj0, ii0, enmb)
# Mean melt pond depth:
hpond1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "hpond_d", jj0, ii0, enmb)
hpond2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "hpond_d", jj0, ii0, enmb)
# Mean pond ice thicknes:
ipond1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "ipond_d", jj0, ii0, enmb)
ipond2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "ipond_d", jj0, ii0, enmb)
# Snow density compaction:
rhoscmp1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "rhos_cmp_d", jj0, ii0, enmb)
rhoscmp2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "rhos_cmp_d", jj0, ii0, enmb)
# Rate of snow loss to leads (kg / (m2*sec))
fsloss1 = read_cicefld_gridpnt(fday1, init_date, init_hr, "fsloss_d", jj0, ii0, enmb)
fsloss2 = read_cicefld_gridpnt(fday2, init_date, init_hr, "fsloss_d", jj0, ii0, enmb)





