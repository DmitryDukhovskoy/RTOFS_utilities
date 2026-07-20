import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import xarray as xr

def find_file(rdate, pthice):
  """
    Find GLORYS file name given rdate = YYYYMMDD and GLORYS ice path
  """
  from pathlib import Path
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


