import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import xarray
from copy import copy
import matplotlib.colors as colors

def sens_tests_info(enmb):
  """
    Sensitivity experiments with datm UFS 
    ai - ice conc, hi - ice thickn, hs - snow thickn
    qi - ice enthalpy adjusted in 1st layer to match surf T / or snow T 
    thermo - adjusted ice conudct --> bubbly and dSdt_slow_mode S relax in ice
    fbrd - adjust snow-ice freeboard to keep snow-ice intrf >= sea level
    ITDrdg - snow distribution, snow phys on
    smtrphs - snow metamorphysm is on

    experiments 0 - 25 - January 2025, 14-day f/casts initialized from GDAS analysis
                         for GFSv17 f/casts

    experiments < 20 - corrections implemented in the S. Ocean only
    epxeriments >= 20 are with both N. and S. poles corrected fields (regn="global") 
                      Note that similar runs in old experiments may not be identical due 
                      improvements implemented at a later stage, e.g. inserted snow in old case
                      will be different to the new simulation because added freeboard and ice 
                      enthalpy adjustments

  experiments 30-35 - July 4, 2025, 28 days forecasts, use hsU for all hsnow (hs=hsU)

  hsU - designates updated hsnow insertion code (see extp 11: ice freeboard + enthalpy)
  hsU = hs+qi+fbrd  - updated hsnow code
  hs0 - original code without ice enth. and ice freeboard correction
  SNPHYS = thermo+ITDrdg+snphys 

  experiments 40 - 42 - init 2024/07/01 using SFS initial conditions
                        to compare with SFS runs

  epxeriment 45-47 - 2025/07/04 using RTOFS ice / snow IC
                     insert hsnow - use updated code (=hsU)
  """
  EXPTS = {
    "01" : "control",                         # control January 2025  
    "02" : "ai+hi+hs0+qi+fbrd+thermo",
    "03" : "ai+hi+hs0+qi+fbrd+thermo+ITDrdg",
    "04" : "hs0",                             # deleted - use 21
    "05" : "ai+hi+hsU+thermo+ITDrdg+snphys",
    "06" : "EMPTY",
    "07" : "ai",
    "08" : "ai+hs0",
    "09" : "ai+hi",
    "10" : "ai+hi+hs0",
    "11" : "ai+hi+hsU",                  # ai+hs+(hs+qi+fbrd)
    "20" : "ai",                         # global: N. and S. poles
    "21" : "hsU",                        # global
    "22" : "ai+hsU",
    "23" : "ai+hi",
    "24" : "ai+hi+hsU",
    "25" : "ai+hi+hsU+SNPHYS",
    "30" : "control",                    # control July 4 2025
    "31" : "ai",                         # iconc 07/04/2025
    "32" : "ai+hi",
    "33" : "ai+hs",
    "34" : "ai+hi+hs",
    "35" : "ai+hi+hs+SNPHYS",
    "40" : "control",                  # July 1 2024 - similar to SFS runs
    "41" : "ai",                       # similar to SFS ai, 07/01/2024
    "42" : "ai+hi",                    # similar to SFS ai+hi 07/01/2024
    "45" : "RTOFS ai+hi",              # RTOFS IC July 4 2025
    "46" : "RTOFS ai+hi+hs",            # RTOFS IC July 4 2025
    "47" : "RTOFS ai+hi+hs+SNPHYS",
   }

  key = f"{enmb:02d}"    
  if key in EXPTS:
    sinfo = EXPTS.get(key) 
  else:
    raise Exception(f"key is not found for {enmb}, add experiment ...")

  return sinfo

def sfs_tests_info(enmb):
  """
    Sensitivity experiments with SFS configuration 
    ai - ice conc, hi - ice thickn, hs - snow thickn
    qi - ice enthalpy adjusted in 1st layer to match surf T / or snow T 
    thermo - adjusted ice conudct --> bubbly and dSdt_slow_mode S relax in ice
    fbrd - adjust snow-ice freeboard to keep snow-ice intrf >= sea level
    ITDrdg - snow distribution, snow phys on
    smtrphs - snow metamorphysm is on

    Tfrz - updated SST in MOM6 = Tfrz * aice 

    control runs use SFS IC (from CPC: iconc, ithkn, etc.)

    2024/07/01 - for varying time period
    expt01 and expt02

  2025/07/01 - all expts with MOM6 Tfrz (except for control run)
  RTOFS ice thickness inserted, NSIDC ice conc, hsnow - clim
  expt > 02
  all experiments use updates in snow enthalpy and ice/snow interface updated (to keep aboce sea level)

  hsU - designates updated hsnow insertion code (see extp 11 in DATM UFS: ice freeboard + enthalpy)
  hsU = hs+qi+fbrd  - updated hsnow code
  ALL = ai+hi+hsU+ITDrdg+snphys

  in all expts with sn. physics: thermo 
  SNPHYS = thermo+ITDrdg+snphys 
  radpar - changed parameters for dlt Eddington to increase snow albedo during melting
  sealvl - ice melt pond sea level parameterization 

  Switching to experiments with standard setup:
  ai+hi+hsU+ITDrdg+snphys+radpar with MOM6 Tfrz = radpar with tr_pond_lvl parameterization

  SFSb2 - SFS beta2 parameters: pond topo, + ALL initializations
  """
  EXPTS = {
    "00" : "Sat.clim",          # reserved field for satellite-derived climatology (ice thkn or snow)
    "01" : "control 20240701",
    "02" : "ai+hi+hsU+ITDrdg+snphys+Tfrz",
    "03" : "control 20250701",    # init ice - from CPC used in SFS
    "04" : "ai+hi",
    "05" : "ALL+radpar+sealvl",      # sealvl pond  with tuned pond param: pndaspect=1.2, apnd_sl=0.2
    "06" : "ALL+radpar+pondtopo",    # expt07 but pond parameterization topo
    "07" : "ALL+radpar+pondlvl",     # pond_lvl with all other settings 
    "08" : "ALL+radpar+pondtopo+snclim",     # pond_topo with updated summer snow clim
    "09" : "ALL+radpar+nopond",
    "10" : "ALL+radpar+sealvl+snclim",        # sea level pond parameterization
    "11" : "ALL+radpar+sealvl+drain+snclim",  # sea level pond + change tscale_pnd_drain = 0.5 (default=10)
    "12" : "e11+MLgbr3-nsidc",    # similar to expt11 with ML ithkn (input: iconc NSIDC + GLORYS)
    "13" : "SFSb2+MLgbr3-nsidc",  # SFSbeta2 params + ALL (similar to expt08),  ML ice thickness (iconc NSIDC + GLORYS)
    "14" : "SFSb2+MLgbr3_20250701", # SFSbeta2 params + ALL,  ML ice thickness input: GDAS, lower R_snw (1.8)
    "61" : "ai+hi+hsU+ITDrdg+snphys_debug",   # cat. ice output for debugging
    "62" : "ai+hi+hsU+ITDrdg+snphys",            #  same as 61, every step output 1 day run
    "63" : "ai+hi+hsU+ITDrdg+snphys+pondpar",   #  tr_pond_lvl params changed, every step out 
    "64" : "ai+hi+hsU+ITDrdg+snphys+pondpar+radpar",      #  tr_pond_lvl params changed + rad param (expt07), every step out
    "65" : "ai+hi+hsU+ITDrdg+snphys+sealvl+radpar",      #  tr_pond_seallvl + all from expt 64
   }

  key = f"{enmb:02d}"    
  if key in EXPTS:
    sinfo = EXPTS.get(key) 
  else:
    raise Exception(f"key is not found for {enmb}, add experiment ...")

  return sinfo


def sens_tests_colors():
  # Line colors:
  CLRS      = np.array([
      [0.00, 0.45, 0.70],  # blue
      [0.90, 0.17, 0.31],  # red
      [0.00, 0.62, 0.38],  # green
      [0.90, 0.60, 0.00],  # orange
      [0.95, 0.90, 0.25],  # yellow
      [0.80, 0.47, 0.65],  # pink/violet
      [0.35, 0.70, 0.90],  # light blue
      [0.50, 0.39, 0.64],  # purple
      [0.55, 0.63, 0.79],  # steel blue
      [0.40, 0.60, 0.00],  # olive green
      [0.70, 0.30, 0.00],  # brown
      [0.70, 0.00, 0.30],  # wine red
      [0.00, 0.55, 0.75],  # teal
      [0.75, 0.75, 0.75],  # light gray
      [0.30, 0.30, 0.30],  # dark gray
      [0.10, 0.80, 0.60],  # aqua green
      [0.55, 0.20, 0.60],  # plum purple
      [0.20, 0.70, 0.30],  # jade green
      [0.80, 0.55, 0.35],  # tan
      [0.25, 0.25, 0.55]   # deep indigo
  ])
  return CLRS

def gfs_retro_runs(run_name, node_nm, model='ice'):
  """
    Return list of available retro runs  
    Saved on current machine
  """
  if model == 'ice':
    if node_nm == 'gaea':
      match run_name:
        case "retrov17_01_stream4":
          RUNS=['2025060406', '2025060412', '2025060418', '2025060500', '2025060506']
        case "retrov17_01_stream1a":
          RUNS=['2022090900']
        case _:
          print(f"unrecognized {run_name}, {node_nm}, or model={model}")
          raise Exception("Check inputs for gfs retro experiments")

  return RUNS


def plot_polar2d(A2d, hlon, hlat, regn='north', fgn=1, clrname='ice_thkn', 
                 rmin=None, rmax=None, sttl='Field A2d',
                 btx="mod_gfs_cice_anls.py"):
  """
    Quick stereorgraphic map of a 2D field (A2d)
    hlon, hlat - geogra coord, model grid

  """
  import mod_colormaps as mclrmps
  from mod_utils_fig import minmax_clrmap, colorbar_horiz, bottom_text
  from mpl_toolkits.basemap import Basemap, cm
  
  if clrname == 'ice_thkn':
    clrmp = mclrmps.colormap_ice_thkn()
  elif clrname == 'ice_conc':
    clrmp = mclrmps.colormap_ice_conc()
    
  if rmin is None or rmax is None:
    rmin, rmax = minmax_clrmap(A2d)

  clrmp.set_bad(color=[0.1, 0.1, 0.1])
  cntr_clr = [0.9,0.,1]
    
  if regn == 'south':
    m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l')
    parallels = np.arange(-80,-10,10.)
    meridians = np.arange(-360,359.,45.)
  elif regn == 'north': 
    m = Basemap(projection='npstere',boundinglat=60,lon_0=-10,resolution='l')
    parallels = np.arange(50, 90, 5)
    meridians = np.arange(-360, 359., 45.)

  xh, yh = m(hlon,hlat) # GFS coords

  plt.ion()
  
  fig1 = plt.figure(fgn, figsize=(9,8))
  plt.clf()        
  ax1 = plt.axes([0.1, 0.12, 0.8, 0.8])
                   
  if regn == 'north':
    m.drawparallels(np.arange(60, 90, 10), labels=[0,0,0,0])
  elif regn == 'south':
    m.drawparallels(np.arange(-80, -50, 10), labels=[0,0,0,0])
                    
  m.drawmeridians(np.arange(-180, 180, 45), labels=[0,0,0,0])
  m.drawcoastlines()
  
  img = ax1.pcolormesh(xh, yh, A2d, cmap=clrmp, vmin=rmin, vmax=rmax)
  
  ax1.set_title(sttl, fontsize=12)

  #Plot colrbar
  clb = colorbar_horiz(fig1, ax1, img, rmin=rmin, rmax=rmax, decim=2, extd='max')

  fig1.canvas.draw()

  pos_clb = clb.ax.get_position()
  bot_clb = pos_clb.y0 
  pbtm = bot_clb - 0.05

  bottom_text(btx, pos=[0.02,0.02], fsz=8)

  return fig1, ax1, img, clb


