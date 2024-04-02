"""
  Find experiment # for GOFS3.1 reanalysis 
  epxeriments 53.X
  for given date
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys


PPTHN = '/home/Dmitry.Dukhovskoy/python'
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')

import mod_time as mtime


def gofs31_expt53X(dnmb):
  nexpt = 53.0

  d30 = mtime.datenum([1999,3,31])
  d31 = mtime.datenum([2000,12,31])
  d32 = mtime.datenum([2003,6,30])
  d33 = mtime.datenum([2005,6,30])
  d34 = mtime.datenum([2007,6,30])
  d35 = mtime.datenum([2009,6,30])
  d36 = mtime.datenum([2011,6,30])
  d37 = mtime.datenum([2013,6,30])
  d38 = mtime.datenum([2014,12,30])
  d39 = mtime.datenum([2015,12,31])
  if dnmb > d30 and dnmb <= d31:
    nexpt = 53.1
  elif dnmb > d31 and dnmb <= d32:
    nexpt = 53.2
  elif dnmb > d32 and dnmb <= d33:
    nexpt = 53.3
  elif dnmb > d33 and dnmb <= d34:
    nexpt = 53.4
  elif dnmb > d34 and dnmb <= d35:
    nexpt = 53.5
  elif dnmb > d35 and dnmb <= d36:
    nexpt = 53.6
  elif dnmb > d36 and dnmb <= d37:
    nexpt = 53.7
  elif dnmb > d37 and dnmb <= d38:
    nexpt = 53.8
  elif dnmb > d38 and dnmb <= d39:
    nexpt = 53.9

  cexpt = f'{nexpt:3.1f}'

  return nexpt, cexpt

