"""
  Read binary file with computed Rossby radius
"""
import os
import numpy as np
import struct
import pickle
import sys
import pdb
import importlib
import matplotlib.pyplot as plt
sys.path.append('/data/home004/dmitry.dukhovskoy/python_codes/MyPython')
sys.path.append('/data/home004/dmitry.dukhovskoy/python_codes/read_wod')
import mod_utils_fig
# importlib.reload(mod_utils_fig)
from mod_utils_fig import bottom_text


YR1 = 2000
YR2 = 2001

f_addWOA = 0
if f_addWOA > 0:
  sfx = '_woa'
else:
  sfx = ""

pthrsb = '/data/ncei2/OCL/sys.inf/'  # dir for Rossby r. output
frsb = pthrsb+'RossbyRadius'+sfx+'.s'
pthout = '/data/ncei2/work/dmitry.dukhovskoy/data_output/'

btx = 'anls_processed_stat.py'


PRBS=["CTD","OSD","MRB","PFL","DRB"]
ALL = {}
for prbnm in PRBS:
  qcnm = prbnm+'QC'
  ALL[qcnm]=[]


for YR in range(YR1,YR2+1):
  fout = pthout+'processed_casts{0}{1}.pkl'.format(YR,sfx)
  print('Loading saved processed records #'+fout)
  with open(fout,'rb') as fid:
    PRCSD = pickle.load(fid)

  keysP = list(PRCSD.keys())

  for prbnm in PRBS:
    qcnm = prbnm+'QC'
    try:
      ix = keysP.index(qcnm)
    except:
      ix = -1
      continue

    dmm = PRCSD.get(qcnm)
    dmm0= ALL.get(qcnm).copy()
    if dmm:
      dmm0.extend(dmm)

    ALL.update({qcnm: dmm0}) 
  

# Hist of QC flags
plt.ion()
if f_addWOA == 1:
  ctl = "WOD Rossby Radius - WOA clim filled deep gaps, km"
  ctl2= "QC flags, used WOA clim" 
  ctl3= "Fraction of QC-passed T/S profiles, used WOA"
else:
  ctl = "WOD Rossby Radius - no WOA clim, km"
  ctl2= "QC flags, no WOA clim" 
  ctl3= "Fraction of QC-passed T/S profiles, no WOA"

# Plot QC flags
plt.close('all')
fig, axs = plt.subplots(3,2,num=1,figsize=(8,9),clear=True, \
           sharex='col', sharey='row')
fig.suptitle(ctl2)
icc = -1
nQC=11     # # of QC flags + 1 for 0=ok
for prbnm in PRBS:
  qcnm = prbnm+'QC'
  qcf = ALL.get(qcnm).copy()
  icc += 1
  ax = axs.flat[icc]
  ax.set_title(prbnm)

  if not qcf:
    continue

  ax.hist(qcf,edgecolor=[0,0,0],range=(0,nQC),bins=nQC,density=True)
  ax.xaxis.set_ticks(np.arange(0,nQC+1,1))
  for ax in axs.flat:
    ax.label_outer()

#  for ax in axs.flat:
#      ax.set(xlabel='x-label', ylabel='y-label'

Txt = ['QC flags:']
Txt.append(' 0 = QC passed')
Txt.append(' 1 = not enough observed values')
Txt.append(' 2 = big difference in # obs. T/S prof')
Txt.append(' 3 = profiles too shallow wrt to local depth')
Txt.append(' 4 = shallow depth')
Txt.append(' 5 = no obs in the upper ocean')
Txt.append(' 6 = top obs. nans, cannot be reconstructed')
Txt.append(' 7 = bottom obs. nans, cannot be reconstructed')
Txt.append(' 8 = last obs too far from bottom')
Txt.append(' 9 = large "jumps" in T/S profiles')
Txt.append('10 = Error in lon/lat from header')
nt = len(Txt)
ax = axs.flat[5]
x0 = 0.0
y0 = 1.0
for itt in range(nt):
  ssl = Txt[itt]
  x1  = x0
  y1  = y0-0.07*(itt+1)
  ax.text(x1,y1,ssl)


bottom_text(btx,fsz=9)
  

# For all selected casts for processing, 
# Plot QC rejected/accepted casts for each probe
Rscs = np.zeros(len(PRBS))
ii = -1
for prbnm in PRBS:
  qcnm = prbnm+'QC'
  qcf = ALL.get(qcnm).copy()
  Qcf = np.array(qcf)
  ntot=len(qcf)
  npass=len(np.where(Qcf==0)[0])

  ii += 1
  if ntot > 0:
    Rscs[ii] = float(npass)/float(ntot)


nprb = len(PRBS)
xx = np.linspace(1,nprb,nprb)

fig2 = plt.figure(2,figsize=(8,6), constrained_layout=False)
fig2.clf()
ax1 = plt.axes([0.1, 0.15, 0.8, 0.7])
bw = 0.8
bclr = [0.6, 0.7, 1.]
ax1.bar(PRBS,Rscs,bw,color=bclr)

ctl3 = '{0} {1}-{2}'.format(ctl3,YR1,YR2)
ax1.set_title(ctl3)

bottom_text(btx,fsz=9)


