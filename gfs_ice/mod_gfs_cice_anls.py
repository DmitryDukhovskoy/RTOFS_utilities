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
  """
  EXPTS = {
    "01" : "control",
    "02" : "ai+hi+hs+qi+fbrd+thermo",
    "03" : "ai+hi+hs+qi+fbrd+thermo+ITDrdg",
    "04" : "hs",
    "05" : "ai+hi+hs+qi+fbrd+thermo+ITDrdg+smtrphs",
    "06" : "  ",
    "07" : "ai",
    "08" : "ai+hs",
    "09" : "ai+hi",
    "10" : "ai+hi+hs",
    "11" : "ai+hi+hs+qi+fbrd",
   }

  key = f"{enmb:02d}"    
  sinfo = EXPTS.get(key) 

  return sinfo


def sens_tests_colors():
  # Line colors:
  CLRS      = np.array([
      [0.00, 0.45, 0.70],  # blue
      [0.90, 0.17, 0.31],  # red
      [0.00, 0.62, 0.38],  # green
      [0.95, 0.90, 0.25],  # yellow
      [0.80, 0.47, 0.65],  # pink/violet
      [0.90, 0.60, 0.00],  # orange
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
