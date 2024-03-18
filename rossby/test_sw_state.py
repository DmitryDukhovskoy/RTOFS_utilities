"""
  Test Sea Water State equations
"""
import os
import numpy as np
import matplotlib.pyplot as plt
#import torch
import sys
import pdb
import mod_swstate
#importlib.reload(mod_swstate)
from mod_swstate import sw_press
from mod_swstate import adiab_Tgrad
from mod_swstate import sw_ptmp
from mod_swstate import sw_dens0
from mod_swstate import sw_smow
from mod_swstate import sw_seck
from mod_swstate import sw_dens

def test_sw(t=20.,s=34.,pr1=1000.,pref=2000.,p=5000.)
  """
    t - water temp
    s - water salinity
    pr1 - pressure at depth 1 (dB)
    pref - pressure reference for computing sigma_pref - density (dB)
    p  - pressure dB for secant bulk formula
  """

  t=20.0
  s=34.0
  pres=10.0
  adtg = adiab_Tgrad(s,t,pres)
  print('adiab T grad(s={0},t={1},p={2}db) = {3}'.
         format(s,t,pres,adtg))

  # Test: potential T from pressure P1 to Pref
  # for t=20, s=34, pr1=1000db, pref=2000db
  # pot. temp = 20.1958294083943 ...
#  t=20.0
#  s=34.0
#  pr1=1000.0
#  pref=2000.0
  ptmp = sw_ptmp(s,t,pr1,pref)
  print('Potential T(s={0},t={1},p1={2}db, pref={3}) = {4}'.
         format(s,t,pr1,pref,ptmp))


# Test rho pure water at P=0
# at t=24.4 rho_S0 = 997.198518670105
#  t=24.4
  rho_S0 = sw_smow(t)
  print('Pure water dens(T={0}) = {1}'.format(t,rho_S0))

  # Test depth ---> pressure conversion
  # should get 7500.00 db
  #pres_db, pres_pa =sw_press(7321.45,30)

#
# Test sea water dens at the surface P=0 
# t=11.5, s=45.2, rho=1034.63853438209
#  t=11.5
#  s=45.2
  rho_P0 = sw_dens0(s,t)
  print('Surface Water dens(S={0}, T={1}) = {2}'.format(s,t,rho_P0))


# Secant bulk formula
# t=11.5, s=45.2, P=5000 db
# K secnat = 25057.5763828814
#  t=11.5
#  s=45.2
#  p=5000.0
  Kscnt = sw_seck(s,t,p)
  print('Secant bulk modulus(S={0}, T={1}, P={2}) = {3} bar'.format(s,t,p,Kscnt))

#
# in situ water dens
# t=11.5, s=45.2, P=5000 db
# 1055.7041012411848 kg/m3
#  t=11.5
#  s=45.2
#  p=5000.0
  rho = sw_dens(s,t,p)
  print('rho(S={0}, T={1}, P={2}) = {3} kg/m3'.format(s,t,p,rho))

  return

