"""
  Read/plot river  rtofs output fields
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.colors import ListedColormap
from copy import copy

#from mpl_toolkits.basemap import Basemap, shiftgrid

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')

rdate   = '20220617'
#pthbin  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_para7b/hycom/rtofs.'\
#           +rdate+'/'
pthbin  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/rtofs_paraCd/fix/'
pthgrid = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/hycom_fix/'

mo_read = 3  # month to read

flnm = 'forcing.rivers'
fina = pthbin+flnm+'.a'
finb = pthbin+flnm+'.b'

ftopo = 'regional.depth'
fgrid = 'regional.grid'


import mod_read_hycom
#importlib.reload(mod_read_hycom)
from mod_read_hycom import read_grid_topo, read_hycom, read_topo
from mod_read_hycom import dx_dy
from mod_utils_fig import bottom_text
import mod_utils as mutil


print('Processing '+fina)

IDM = 4500
JDM = 3298

LON, LAT, HH = read_grid_topo(pthgrid,ftopo,fgrid)

huge = 2.**90
IJDM = IDM*JDM
npad = 4096-IJDM%4096
# Read river:
try:
  fgb = open(fina,'r')
except:
  print('Could not open '+fina)

fga = open(fina,'rb')
F = []

k0 = mo_read-1
fga.seek(k0*(npad+IJDM)*4,0)
dmm = np.fromfile(fga, dtype='>f',count=IJDM) # read 1 layer
amin = np.min(dmm[np.where(dmm<huge)])
amax = np.max(dmm[np.where(dmm<huge)])
print('Reading rivers month={0} min={1:12.8f} max={2:12.8f} m/sec'.\
        format(mo_read,amin,amax))
Friv = dmm.reshape((JDM,IDM))
fga.close()


# Convert to m3/mo
DX, DY = dx_dy(LON,LAT)
Acell = DX*DY
Friv  = Friv*Acell  # m3/sec 



LMSK = HH.copy()
LMSK = np.where(LMSK<0.,0.,1.)
lcmp = mutil.clrmp_lmask(2)

il1 = 2380
il2 = 2480
jl1 = 1860
jl2 = 1910

Friv = np.where(HH>0.,-9.e20,Friv)

# =================
#  PLOT RIVERS Mississippi
# =================
plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.2, 0.8, 0.7])
#ax1.pcolormesh(LMSK, shading='flat',\
#                cmap=lcmp, \
#                vmin=0, \
#                vmax=1)

clrmp = copy(plt.cm.cubehelix_r)
clrmp.set_bad(color=[0.8,0.8,0.8])
clrmp.set_under(color=[0.8,0.8,0.8])
im1 = ax1.pcolormesh(Friv, shading='flat',\
               vmin=0, vmax=1000.,cmap=clrmp)

ax1.set_xlim([il1,il2])
ax1.set_ylim([jl1,jl2])

# Total runoff Mississip:
daa = Friv[1870:1900,2420:2467]
daa[np.where(daa<0.)] = np.nan
Rtot = np.nansum(daa)


ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
             ax1.get_position().y0,0.02,
             ax1.get_position().height])
clb = plt.colorbar(im1, cax=ax2, extend='max')
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
clb.ax.set_yticklabels(ticklabs,fontsize=10)

plt.sca(ax1)
stl = ('River runoff in RTOFS/HYCOM grid, m3/sec, Mo={0}, Tot={1:7.1f}'.\
        format(mo_read,Rtot))
ax1.set_title(stl)

btx = 'plot_rivers.py'
bottom_text(btx)


