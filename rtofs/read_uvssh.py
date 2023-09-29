"""
  Example: Read staggered u, v fields from RTOFS *archv* binary
           Read ssh
"""
import os
import numpy as np
import sys
import importlib

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

import mod_read_hycom as mhycom


rdate  = '20230129'
pthbin = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/rtofs.'+\
          rdate+'/ocnqc_logs/profile_qc/'
pthhcm = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/rtofs.'+rdate+'/'
flhcm  = 'rtofs_glo.t00z.n-24.archv'

fina = pthhcm+flhcm+'.a'
finb = pthhcm+flhcm+'.b'

huge = 1.e20
rg   = 9806.
IDM, JDM, KDM = mhycom.hycom_dim(fina,finb)


# Read U-barotropic:
F,nn,mm,ll = mhycom.read_hycom(fina,finb,'u_btrop')
F[np.where(F>huge)] = np.nan
Ub = F.copy()

# Read baroclinic U component and add barotropic
U   = []
fld = 'u-vel.'
for kk in range (1,KDM+1):
  F,nn,mm,ll = mhycom.read_hycom(fina,finb,fld,rLayer=kk)
  F[np.where(F>huge)] = np.nan
  F = F+Ub
  if kk == 1:
    U = np.copy(F)
  else:
    U = np.dstack((U,F))


# Read ssh:
SSH,nn,mm,ll = mhycom.read_hycom(fina,finb,'srfhgt')
SSH[SSH>huge] = np.nan
SSH = SSH/9.806



