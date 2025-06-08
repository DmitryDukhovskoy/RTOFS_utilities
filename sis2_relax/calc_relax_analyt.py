"""
  You start with
  (X_new - X_old)/dt = rlx*(X_an - X_new)

  you can massage it into:
  X_new (1 + dt*rlx) = X_old + dt*rlx*X_an

  then:
  X_new = 1/(1 + dt*rlx) * (X_old + dt*rlx*X_an)
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


Trlx_hr = 0.08
Trlx_sec = Trlx_hr*3600
Irlx = 1/Trlx_sec

dt_slow = 3600. # slow time step, sec
dt = dt_slow
#Iresttime = 0.4557e-4
Iresttime = Irlx
damp = dt * Iresttime
I1pdamp = 1./(1. + damp)
Xold = 0.9796
Xref = 0.8841

Xnew = I1pdamp * (Xold + Xref * damp)

print(f'Xnew = {Xnew}')

