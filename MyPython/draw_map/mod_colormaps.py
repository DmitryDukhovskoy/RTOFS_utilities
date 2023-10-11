"""
  Colormaps
"""
import numpy as np
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.colors import ListedColormap


def clrmp_BlRd(Ncmp):
# Blue - white - Red:
  CLR =[[14.2,        0.0,       85.0,    1],
            [   9.7,       11.2,       97.0,    1],
            [ 16.0,       34.2,      100.0,    1],
            [ 24.0,       53.1,      100.0,    1],
            [ 34.0,       69.2,      100.0,    1],
            [ 46.0,       82.9,      100.0,    1],
            [ 60.0,       92.0,      100.0,    1],
            [ 74.0,       97.8,      100.0,    1],
            [ 92.0,      100.0,      100.0,    1],
            [100.0,      100.0,       92.0,    1],
            [100.0,       94.8,       74.0,    1],
            [100.0,       84.0,       60.0,    1],
            [100.0,       67.6,       46.0,    1],
            [100.0,       47.2,       34.0,    1],
            [100.0,       24.0,       24.0,    1],
            [ 97.0,       15.5,       21.0,    1],
            [ 85.0,        8.5,       18.7,    1],
            [ 65.0,        0.0,       13.0,    1]];

  CLR = np.array(CLR)
  CLR = CLR/100.
  CLR[:,3] = 1.
  CMP = create_colormap(CLR,Ncmp)
  return CMP


def clrmp_BlGrRd(Ncmp):
# Blue - lightGreen - Red:
  CLR =[[14.2,        0.0,       85.0,    1],
            [   9.7,       11.2,       97.0,    1],
            [ 16.0,       34.2,      100.0,    1],
            [ 24.0,       53.1,      100.0,    1],
            [ 34.0,       69.2,      100.0,    1],
            [ 46.0,       82.9,      100.0,    1],
            [ 60.0,       92.0,      100.0,    1],
            [ 74.0,       97.8,      100.0,    1],
            [ 74.0,      100.0,       92.0,    1],
            [ 92.0,      100.0,       74.0,    1],
            [100.0,       94.8,       74.0,    1],
            [100.0,       84.0,       60.0,    1],
            [100.0,       67.6,       46.0,    1],
            [100.0,       47.2,       34.0,    1],
            [100.0,       24.0,       24.0,    1],
            [ 97.0,       15.5,       21.0,    1],
            [ 85.0,        8.5,       18.7,    1],
            [ 65.0,        0.0,       13.0,    1]];

  CLR = np.array(CLR)
  CLR = CLR/100.
  CLR[:,3] = 1.
  CMP = create_colormap(CLR,Ncmp)
  return CMP


def clrmp_Nvalues(Ncat,Ncmp):
# Colormap for distinct groups
# Specify N main colors
# Specify Ncmp - total # of colorshades

  CLR0 =[[0.95,  0.95, 1],  # cat 1 
         [ 0.,  0.,  0.5],
         [0.95, 1, 0.95],
         [0.,  0.5,  0],    # cat 2
         [1., 0.95, 0.95],
         [0.5, 0., 0.],     # cat 3
         [1., 0.95, 1.],
         [0.5, 0., 0.5],    # cat 4
         [0.95, 1., 1.],
         [0., 0.5, 0.5],    # cat 5
         [1., 1., 0.95],
         [0.5, 0.5, 0.],    # cat 6
         [1., 0.95, 0.9],
         [0.75, 0.5, 0],    # cat 7
         [0.9, 1., 0.95],
         [0., 0.75, 0.5],  # cat 8
         [0.95, 1., 0.9],
         [0.5, 0.75, 0.],
         [1., 0.92, 0.95],
         [0.7, 0.3, 0.5]]

  CLR0 = np.array(CLR0)
  ny, nx = CLR0.shape
  a = np.insert(CLR0,nx,1.0, axis=1)

  if Ncat > ny:
    print("clrmp_Nvalues: WARNING # specified categ {0} > # of color groups {1}".format(Ncat,ny))

  if Ncat < ny:
    CLR = CLR0[:2*Ncat]
  else:
    CLR = CLR0.copy()

  CMP = create_colormap(CLR, Ncmp)
  return CMP



def create_colormap(CLR, Ncmp, cmp_obj=True):
# Mix main colors in CLR adding shadings
# Transitioning from CLR(i) to CLR(i+1)
# Ncmp - total # of colors in the colormap
# In most cases Ncmp > length(CLR)
# otherwsie CLR is returned with no changes
# cmp_obj - True ==> returns as colormap object
#           False --> returns as RGB array
  nClr = CLR.shape[0]

# If number of shades is less or eq. the # of Main colors 
# use the Colors and no shades
  if Ncmp < nClr:
    print('create_colormap:')
    print('Specified N of colors {0} </= N of Main Colors {1}'.format(Ncmp,nClr))
    print(' Adjust to Ncmp clrs to {0}'.format(nClr))
    print(' Colormap not changed')
    vals = CLR
    CMP = ListedColormap(vals)
    return CMP

#
# Define # of colors for each color shade
# Then linearly interpolate btw the colors
  nInt = nClr-1
  newCLR = np.empty((0,4))
  for ii in range(nInt):
    clr1 = CLR[ii,:]
    clr2 = CLR[ii+1,:]
# Last interval - no extra shade
    if ii < nInt-1:
      nShades = np.int(np.round(Ncmp/nInt))+1
    else:
      nShades = np.int(np.round(Ncmp/nInt))
    vals = np.ones((nShades,4))
    for jj in range(3):
      vals[:,jj] = np.linspace(clr1[jj],clr2[jj],nShades)

# Delet last row - same as the next color
    nv = vals.shape[0]
    if ii < nInt-1:
      vals = vals[0:nv-1]
    newCLR = np.append(newCLR,vals, axis=0)

  if cmp_obj:
    CMP = ListedColormap(newCLR)
  else:
    CMP = np.array(newCLR)
# breakpoint()

  return CMP

def colormap_conc():
  """
   Prepare colormap for sea ice conc maps
   Creates ListedColormap object
   to get colors:
   cmpice.colors
   cmpice.N - # of colors
  
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
         [8,    80, 174],
         [255, 255, 255]]
  
  CLR = np.array(CLR)/255.
  CLR = np.flip(CLR, axis=0)
  CMP = mclrs.create_colormap(CLR, 200)
  
  return CMP

