"""
  Plot interpolated TOPO ETOPO2 ---> WOA0.25 degree
  Using Fortran code: /data/home004/dmitry.dukhovskoy/fortran etopo2woa025.F90
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('/data/home004/dmitry.dukhovskoy/python_codes/MyPython')
import mod_utils_fig
from mod_utils_fig import bottom_text
from mod_utils_fig import plot_fld2D
from copy import copy
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

pthout = '/data/ncei2/work/dmitry.dukhovskoy/topo/woa/'
fout   = pthout+'etopo2woa25_invdst.dat'

# Read topo:
# Reading fortran binary:
#
#with open(fout,'rb') as fid:
fid = open(fout,'rb')
try:
  fid.seek(0)
  imm = np.fromfile(fid,dtype='<i4', count=1) # bytes in following record
  idm = np.fromfile(fid,dtype='<i4', count=1)[0]
  imm = np.fromfile(fid,dtype='<i4', count=1) # bytes in following record

  imm = np.fromfile(fid,dtype='<i4', count=1) 
  jdm = np.fromfile(fid,dtype='<i4', count=1)[0]
  imm = np.fromfile(fid,dtype='<i4', count=1) 

  ijdm = idm*jdm
  imm = np.fromfile(fid,dtype='<i4', count=1) 
  lonW = np.fromfile(fid,dtype='<f4', count=idm)
  imm = np.fromfile(fid,dtype='<i4', count=1) 
 
  imm = np.fromfile(fid,dtype='<i4', count=1) 
  latW = np.fromfile(fid,dtype='<f4', count=jdm)
  imm = np.fromfile(fid,dtype='<i4', count=1) 

  imm = np.fromfile(fid,dtype='<i4', count=1)
  HH  = np.fromfile(fid,dtype='<f4', count=ijdm) 

finally:
 fid.close()

#A0 = np.copy(HH)
HH = HH.reshape(jdm,idm,order='C')


clrmp = copy(plt.cm.jet)
clrmp.set_over(color=[0.5, 0.5, 0.5])

plt.ion()
fig1 = plt.figure(2, figsize=(8,8), constrained_layout=False)
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
cl1=-8000.
cl2=0.

f_indx=1
if f_indx>0:
		im1 = ax1.pcolormesh(HH,shading='flat',cmap=clrmp,
														norm=colors.Normalize(vmin=cl1, vmax=cl2))
		ax1.axis('scaled')
		ax1.set(xlim=(0,idm),ylim=(0,jdm))
else:
		im1 = plt.pcolormesh(lonW,latW,HH,shading='nearest',cmap=clrmp,
															norm=colors.Normalize(vmin=cl1, vmax=cl2))
		ax1.axis('scaled')
		ax1.set(xlim=(-180,180),ylim=(-90,90))

im1.set_clim(cl1,cl2)
ax1.set_title('TOPO WOA0.25 from ETOPO2, inv.dist interp')
ax1.autoscale(enable=True, axis='both', tight=True)

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom",size="5%",pad=0.25)
hb = fig1.colorbar(im1,ax=ax1,orientation='horizontal',cax=cax)

btx = 'plot_topoWOA_Fort.py' 
bottom_text(btx,pos=[0.1,0.1])
#  plt.colorbar(im1)
plt.show()


#
# Plot Topo difference btw 2 methods Inv. Dist and bilinear
fout = pthout+'etopo2woa025_bilinear.dat'
fid = open(fout,'rb')
try:
  fid.seek(0)
  HB = np.fromfile(fid,dtype='<f4', count=ijdm)  

finally:
  fid.close()

HB = HB.reshape(jdm,idm)

dH = HH-HB

clrmp2 = copy(plt.cm.RdBu)
fig2 = plt.figure(3, figsize=(8,8), constrained_layout=False)
plt.clf()
ax1 = plt.axes([0.1, 0.4, 0.8, 0.5])
cl1=-20.
cl2=20.

im2 = ax1.pcolormesh(dH,shading='flat',cmap=clrmp2,
												norm=colors.Normalize(vmin=cl1, vmax=cl2))
ax1.axis('scaled')
ax1.set(xlim=(0,idm),ylim=(0,jdm))

ax1.set_title('Difference in TOPO WOA0.25 from  inv.dist & bilinear interp')
ax1.autoscale(enable=True, axis='both', tight=True)

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom",size="5%",pad=0.25)
fig2.colorbar(im2,ax=ax1,orientation='horizontal',cax=cax,extend='both')

ax2 = plt.axes([0.1, 0.1, 0.4, 0.25])
dH1d = dH.flatten()
ax2.hist(dH1d,edgecolor=[0,0,0],range=(-20,20),bins=40,density=True)
ax2.set_xticks(range(-20,20,5))
ax2.set_title('TOPO difference (Ind.Dist - Bilinear), m')

btx = 'plot_topoWOA_Fort.py'
bottom_text(btx,pos=[0.1,0.03])
#  plt.colorbar(im1)
plt.show()




