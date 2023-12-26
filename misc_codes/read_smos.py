"""
  Read SMOM satellite salinity data
"""
import os
import numpy as np

from netCDF4 import Dataset as ncFile

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/ncoda_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')

from mod_utils_fig import bottom_text


YR   = 2023
MM   = 9
DM   = 8
sfx1 = '075750'
sfx2 = '085110'
dstr = '{0}{1:02d}{2:02d}'.format(YR,MM,DM)
pth  = '/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data/SMOS/' + dstr + \
       '/wtxtbul/satSSS/SMOS/'
flnm = 'SM_OPER_MIR_OSUDP2_{0}T{1}_{0}T{2}_700_001_1.nc'.format(dstr,sfx1,sfx2)
flnm = 'SM_OPER_MIR_OSUDP2_20230908T102802_20230908T112116_700_001_1.nc'
flinp= pth + flnm

print('dir: ' + pth)
print('Opening ' + flnm)

ncdata = ncFile(flinp, 'r')
LON    = ncdata['Longitude'][:].data.squeeze()
LAT    = ncdata['Latitude'][:].data.squeeze()
AA     = ncdata['SSS_corr'][:].data.squeeze()
AA     = np.where(AA<-1., np.nan, AA)

plt.ion()
fig1 = plt.figure(nfg,figsize=(9,9))
plt.clf()
plt.plot(LAT,AA)
plt.title(flnm)


