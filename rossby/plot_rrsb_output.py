"""
  Read binary file with computed Rossby radius
"""
import os
import numpy as np
import struct
import matplotlib.pyplot as plt

f_addWOA = 0
if f_addWOA > 0:
  sfx = '_woa'
else:
  sfx = ""

pthrsb = '/data/ncei2/OCL/sys.inf/'  # dir for Rossby r. output
frsb = pthrsb+'RossbyRadius'+sfx+'.s'

btx = 'plot_rrsb_output.py'

icc  = 0
itot = 0
Rrsb = []
Unmb = []
rec_size = struct.calcsize('f')

"""
fid = open(frsb,'rb')
fid.seek(0,2)
fend = fid.tell()
nrec_tot = int(fend/rec_size)
fid.seek(0)
print('Total # of records={0}'.format(nrec_tot))
try: 
  bd = fid.read(rec_size)
  val = struct.unpack('f',bd)[0]
  if val > 1.e-10:
    icc += 1
    Rrsb.append(val)

    cpos = fid.tell()
    uniqnmb = cpos/rec_size+1
    Unmb.append(uniqnmb)     
    print('{0} U#{1} Rrsb{2:.1f}'.format(icc,uniqnmb,Rrsb))
except:
  print('Error reading file')
finally:
  fid.close()
"""

with open(frsb,'rb') as fid:
# Determine the size of file
  fid.seek(0,2)
  fend = fid.tell()
  nrec_tot = int(fend/rec_size)
  fid.seek(0)

  for itot in range(nrec_tot):
    if itot % 500000 == 0:
      print(' ... processed {0:.2f}%'.format(itot/nrec_tot*100.))

    try:
      bd = fid.read(rec_size)
      val = struct.unpack('f',bd)[0]
      if val > 1.e-6:
        icc += 1
        Rrsb.append(val)

        cpos = fid.tell()
        uniqnmb = cpos/rec_size+1
        Unmb.append(uniqnmb)     
        print('{0} U#{1} Rrsb{2:.1f}'.format(icc,uniqnmb,val))
 
    except:
      print('Error reading file')
      break


Rrsb = np.array(Rrsb)

#
# Hist of Rossby R
plt.ion()
plt.close('all')
fig1 = plt.figure(1,figsize=(8,5), constrained_layout=False)
fig1.clf()
ax1 = plt.axes([0.1, 0.15, 0.8, 0.7])
ax1.hist(Rrsb,edgecolor=[0,0,0],range=(0,300),bins=30,density=True)
ax1.set_xticks(np.arange(0,300,step=20))
if f_addWOA == 1:
  ctl = "WOD Rossby Radius - WOA clim filled deep gaps, km"
else:
  ctl = "WOD Rossby Radius - no WOA clim, km"

ax1.set_title(ctl)

bottom_text(btx,fsz=9)




