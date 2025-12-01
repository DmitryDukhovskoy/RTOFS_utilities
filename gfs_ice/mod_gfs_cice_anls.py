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
  """
  EXPTS = {
    "01" : "control",
    "02" : "3x snowf",
    "03" : "10x snowf",
    "04" : "hs",
    "05" : "ai+3x snowf",
    "06" : "ai+10x snowf",
    "07" : "ai",
    "08" : "ai+hs",
    "09" : "ai+hi",
    "10" : "ai+hi+hs",
    "11" : "ai+hi+hs+qice+freeb",
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


    
