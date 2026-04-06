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

  experiments 30-35 - July 2025, 28 days forecasts

  hsU - designates updated hsnow insertion code (see extp 11: ice freeboard + enthalpy)
  hsU = hs+qi+fbrd  - updated hsnow code
  hs0 - original code without ice enth. and ice freeboard correction
  SNPHYS = thermo+ITDrdg+snphys 

  """
  EXPTS = {
    "01" : "control",                         # control January 2025  
    "02" : "ai+hi+hs0+qi+fbrd+thermo",
    "03" : "ai+hi+hs0+qi+fbrd+thermo+ITDrdg",
    "04" : "hs0",
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
    "30" : "control",                    # control July 2025
    "31" : "ai",                         # July, iconc 
    "32" : "ai+hi",
    "33" : "ai+hsU",
    "34" : "ai+hsU+hi",
    "35" : "ai+hi+hsU+SNPHYS",
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
