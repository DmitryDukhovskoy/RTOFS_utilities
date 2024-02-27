"""
  Functions for readin/plotting CICE output
  on 0.08 Global tripolar grid used in RTOFS
"""
import numpy as np

def sub_region2D(AA, region='Arctic'):
  """
    Subsample polar regions from the global grid
    input: AA - 2D scalar field

  """
  imd = AA.shape[1]
  jmd = AA.shape[0]

  if region == 'Arctic':
    jS = 2228
    jE = jmd
  elif region == 'Antarctic':
    jS = 0
    jE = 750

  AS = AA[jS:jE,:]

  return AS 


def colormap_conc():
  """
   Prepare colormap for sea ice conc maps
  """
  import mod_colormaps as mclrs
  CLR = [[238, 226, 215],
         [223, 206, 189],
         [216, 162, 107],
         [208, 131,  54],
         [232, 177,  59],
         [232, 208,  50],
         [250, 241, 110],
         [219, 240,  94],
         [210, 250, 162],
         [157, 246,  97],
         [97,  246, 127],
         [35,  202, 157],
         [122, 238, 246],
         [4,   173, 185],
         [25,  154, 253],
         [8,    80, 174]]

  CLR = np.array(CLR)/255.
  CMP = mclrs.create_colormap(CLR, 200)

  return CMP


def read_iceconc_cdr(fpice, fnmb, xP, yP, lon0=-45., regn='Arctic', landnan=True):
  """
    Data downloaded from NESDIS Polar Watch portal
    https://polarwatch.noaa.gov/catalog/
    This is 25-km sea ice concentration data for the Arctic 
    from the Near Real-Time NOAA/NSIDC Sea Ice Concentration 
    Climate Data Record (G10016_v2). Daily and monthly versions 
    are available with data beginning Jan 2021.

    Assuming both projections are polar azimuthal
    collocate poles: xP, yP - pole coord in main figure

    Data projection: central lat = 90, lon = -45
    Need to rotate to the orientation of the lain figure
    lon0 - central longitude in the main figure

    Inofmration about NSIDC polar projection grids:
    https://nsidc.org/data/user-resources/help-center/guide-nsidcs-polar-stereographic-projection
    
  """
  from netCDF4 import Dataset as NetCDFFile
  import mod_time as mtime

  nci = NetCDFFile(fpice)
  DVf   = mtime.datevec(fnmb)
  
# Convert NOAA NESDIS time to calendar/python
# NOAA - # of sec since01/01/1970 00:00
  TS    = nci.variables['time'][:].data
  Ndays = TS/86400.
  TM    = mtime.datenum([1970,1,1]) + Ndays 
  DVobs = mtime.datevec2D(TM)

# Ignore hours, daily data
  dd = np.sqrt( (DVobs[:,0] - DVf[0])**2 + \
                (DVobs[:,1] - DVf[1])**2 + \
                (DVobs[:,2] - DVf[2])**2 )

  Iday = np.where(dd == 0)[0]
  if len(Iday) == 0:
    raise Exception('Requested date {0}/{1}/{2} not found in ice conc {3}'.\
                    format(DVf[0:2],fpice))

  Xobs = nci.variables['xgrid'][:].data
  Yobs = nci.variables['ygrid'][:].data
  Cice = nci.variables['cdr_seaice_conc'][Iday,:,:].squeeze().data

# Find the NPole index:
# assuming (x=0,y=0) = North Pole
#  import mod_interp1D as mintrp1
#  nx = Xobs.shape[0]
#  ny = Yobs.shape[0]
#  IX = np.arange(nx)
#  IY = np.arange(ny)
#  iPobs = mintrp1.pcws_lagr1(Xobs,IX,0.0)
#  jPobs = mintrp1.pcws_lagr1(Yobs,IY,0.0)

  if landnan:
    Cice[np.where(Cice > 1.0)] = np.nan

  lon0_obs = 45.
  if regn == 'Antarctic':
    lon0_obs = 0.

  dPhi = lon0 + lon0_obs
  if abs(dPhi) < 0.01:
#    return Xobs, Yobs, Cice
    Xr, Yr = np.meshgrid(Xobs,Yobs, indexing='xy')

  else:
    theta = np.radians(dPhi)
    cs, sn = np.cos(theta), np.sin(theta)
    Q = np.array([[cs,  sn],
                  [-sn, cs] ])

    nx, ny = Xobs.shape[0], Yobs.shape[0]
    Xr = np.zeros((ny,nx))
    Yr = np.zeros((ny,nx))
    VV = np.zeros((2,1))
    for ix in range(nx):
      for iy in range(ny): 
        VV = np.array([[Xobs[ix]],
                      [Yobs[iy]]])

        VR = np.dot(Q,VV)

        Xr[iy,ix] = VR[0]
        Yr[iy,ix] = VR[1]

# Offset n pole:
# In NCIDS Polar Stereographic, EPSG:3413
# N. Pole is at X=0, Y=0
#  mx, nx = Xr.shape[0], Xr.shape[1]
#  xmin = np.min(abs(Xr))
#  ymin = np.min(abs(Yr))
#  iP, jP = np.where( (np.abs(Xr) == xmin) & (np.abs(Yr) == ymin) )

#  xP0  = Xr[jP,iP]
#  yP0  = Yr[jP,iP]
#  xoff = xP - xP0
#  yoff = yP - yP0

  Xr = Xr + xP
  Yr = Yr + yP

  return Xr, Yr, Cice








