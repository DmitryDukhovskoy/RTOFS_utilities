import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import xarray as xr
from pathlib import Path

def find_file_v0(rdate, pthice):
  """
    Find GLORYS file name given rdate = YYYYMMDD and GLORYS ice path
  """
  try:
    dflice = next(
        Path(pthice).glob(
            f"*_mean_{rdate}_R*.nc"
        )
    )
    print(f"Found file: {dflice}")
    return dflice
  except StopIteration:
    print(f"No file found for {rdate} in {pthice}")
    return None

def find_file(rdate, pthice):
  """
  Find GLORYS file name given rdate = YYYYMMDD and GLORYS ice path.
  """
  pattern = f"*_mean_{rdate}_R*.nc"
  dflice = next(Path(pthice).glob(pattern), None)

  if dflice is None:
    print(f"No file found for {rdate} in {pthice}")
  else:
    print(f"Found file: {dflice}")

  return dflice

