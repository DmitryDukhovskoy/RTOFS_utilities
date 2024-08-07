"""
  To run specific froecast experiment need to 
  specify year/month start 
  the script will prepare a wrapper xml that will call a template xml
"""
NOT FINISHED - use bash script instead

import datetime as dt
import numpy as np
from pathlib import Path
import os
import importlib
import sys
from yaml import safe_load
import textwrap
import argparse

PPTHN = []
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')
sys.path.append('./seasonal-workflow')
import mod_time as mtime


