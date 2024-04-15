"""
  Colormaps
"""
import numpy as np
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.colors import ListedColormap


def minmax_clrmap(dmm,pmin=10,pmax=90,cpnt=0.01,fsym=False):
  """
  Find min/max limits for colormap 
  discarding pmin and 1-pmax min/max values
  cpnt - decimals to leave
  """
  dmm = dmm[~np.isnan(dmm)]
  a1  = np.percentile(dmm,pmin)
  a2  = np.percentile(dmm,pmax)
  cff = 1./cpnt 
  rmin = cpnt*(int(a1*cff))
  rmax = cpnt*(int(a2*cff))
      
  if fsym and (rmin<0. and rmax>0.) :
    dmm = max([abs(rmin),abs(rmax)])
    rmin = -dmm 
    rmax = dmm
    
  return rmin,rmax

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
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap
  nClr = CLR.shape[0]

# If number of shades is less or eq. the # of Main colors 
# use the Colors and no shades
  if Ncmp <= nClr:
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
      nShades = int(np.round(Ncmp/nInt))+1
    else:
      nShades = int(np.round(Ncmp/nInt))
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
#  import mod_colormaps as mclrs
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
  CMP = create_colormap(CLR, 200)
  
  return CMP

def colormap_salin(nclrs=200, clr_ramp=[1,1,1]):
  """
    Colormap for salinity
    low S value ramp to clr_ramp
  """
  from matplotlib import cm
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap
#  import mod_colormaps as mclrs

  btm = cm.get_cmap('rainbow',nclrs)
  ixtop  = round(nclrs*0.1)-1
  clrbtm = btm(range(nclrs))
  chbtm  = np.zeros((ixtop,4))
#
# Add ramp colors at the bottom of clrbar
#  if add_btm == True:
# Add white at the beginning:
  cxbtm  = clrbtm[0,:]

  chbtm[:,3] = cxbtm[3]

  for ik in range(3):
    cc0 = clr_ramp[ik]
    chbtm[:,ik]  = np.linspace(cxbtm[ik],cc0,ixtop)

  chbtm = np.flip(chbtm, axis=0)
  clrbtm = np.insert(clrbtm,0,chbtm, axis=0)

# Add extra colors at the top for better representation of 
# high-S range
  CLR = [[204,   0,   0],
         [153,   0,   0],
         [153,  76,   0],
         [204, 102,   0],
         [255, 229, 192]]
  CLR = np.array(CLR)/255.
  CLR[np.where(CLR > 1.)] =  1.
  CMP = create_colormap(CLR, ixtop, cmp_obj=False)
  clr_high = CMP[0,:]

  nclrs  = clrbtm.shape[0]
  clrtop = clrbtm[-1,:]
  chtop  = np.zeros((ixtop,4))
  chtop[:,3] = cxbtm[3]
  for ik in range(3):
    cc0 = clr_high[ik]
    chtop[:,ik] = np.linspace(clrtop[ik],cc0,ixtop)

# COmbine high S colors at the end of colormap
  clrbtm = np.append(clrbtm, chtop, axis=0)
  clrbtm = np.append(clrbtm, CMP, axis=0)

  newclrs = clrbtm
  newcmp  = ListedColormap(newclrs)

  return newcmp

def colormap_haline(nclrs=200):
  """
    Colormap for salinity
  """
#  import mod_colormaps as mclrs
  CLR = [[200, 230, 245],
         [33,  17,  128],
         [69,   52, 166],
         [39,   10, 204],
         [20,  101, 191],
         [14,  139, 252],
         [45,  171, 201],
         [45,  181, 170],
         [25,  207, 155],
         [89,  194, 120],
         [5,   179,  56],
         [101, 235, 140],
         [182, 227,  92],
         [175, 217,  24],
         [205, 209,  12],
         [240, 214,  20],
         [240, 182,  20],
         [224, 103,  10],
         [245, 141,  99],
         [255, 236, 236]]

  CLR = np.array(CLR)/255.
  CMP = create_colormap(CLR, nclrs)

  return CMP

def colormap_salin2(nclrs=200):
  """
    Colormap for salinity
    bright color palett
    TODO: improve or delete this colormap
  """
#  import mod_colormaps as mclrs
  CLR = [[255, 255, 255],
         [255, 153, 204],
         [255, 102, 255],
         [255,   0, 255],
         [127,   0, 255],
         [178, 102, 255],
         [204, 153, 255],
         [153, 153, 255],
         [102, 102, 255],
         [0,   0,   255],
         [0,   102, 102],
         [51,  153, 255],
         [153, 204, 255],
         [153, 255, 204],
         [51,  255, 153],
         [0,   153,  76],
         [0,   204,   0],
         [102, 255, 102],
         [155, 255, 155],
         [255, 255, 153],
         [255, 255,   0],
         [204, 204,   0],
         [255, 128,   0],
         [255, 178, 102],
         [255, 204, 153],
         [255, 255, 204],
         [255, 153, 153],
         [255,  51,  51]]

  CLR = np.array(CLR)/255.
  CMP = create_colormap(CLR, nclrs)

  return CMP

def colormap_temp2(nclrs=200):
  """
    Temp colormap based on jet with added purple and light red 
    shades to expand the colormap for better T range
  """
  from matplotlib import cm
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap

  CLR = [[242, 186, 238],
         [184,  48, 173],
         [100,  91, 207],
         [ 19,   6, 186],
         [  0,   0, 190],
         [  0,   0, 250],
         [  0,  31, 255],
         [  0,  82, 255],
         [  0, 133, 255],
         [  0, 185, 255],
         [  5, 236, 241],
         [ 47, 255, 200],
         [ 88, 255, 159],
         [129, 255, 117],
         [170, 255,  76],
         [212, 255,  35],
         [254, 237,   0],
         [255, 190,   0],
         [255, 142,   0],
         [255,  95,   0],
         [255,  47,   0],
         [232,   1,   0],
         [174,   0,   0],
         [222, 170, 164]]

  CLR = np.array(CLR)/255.
  CMP = create_colormap(CLR, nclrs)

  return CMP         

def colormap_temp_coldhot(nclrs=200):
  """
    Temp colormap based on jet with added dark purple purple 
    and light red 
    shades to expand the colormap for better T range for cold/warm waters
  """
  from matplotlib import cm
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap

  CLR = [[ 87,  67, 115],
         [ 57,   5, 130],
         [ 89,  27, 176],
         [117,  26, 247],
         [123,  34, 143],
         [184,  48, 173],
         [100,  91, 207],
         [ 19,   6, 186],
         [  0,   0, 190],
         [  0,   0, 250],
         [  0,  31, 255],
         [  0,  82, 255],
         [  0, 133, 255],
         [  0, 185, 255],
         [  5, 236, 241],
         [ 47, 255, 200],
         [ 88, 255, 159],
         [129, 255, 117],
         [170, 255,  76],
         [212, 255,  35],
         [254, 237,   0],
         [255, 190,   0],
         [255, 142,   0],
         [255,  95,   0],
         [255,  47,   0],
         [232,   1,   0],
         [174,   0,   0],
         [150,  47,  47],
         [173, 112, 112],
         [227, 205, 205]]

  CLR = np.array(CLR)/255.
  CMP = create_colormap(CLR, nclrs)

  return CMP

def colormap_landmask(clr0=[0.,0.,0.], clr1=[1.,1.,1]):
  """
    Colormap for land masks 0=land and 1= ocean
    Default: land = black, ocean = white
  """
  CLR = np.array([clr0, clr1])
  CMP = ListedColormap(CLR)

  return CMP

def colormap_speed(nclrs=200):
  """
    Spped colormap from light blue with all colors to red/brown/ light brown
  """
  from matplotlib import cm
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap

  CLR = [[237, 239, 242],
         [206, 214, 237],
         [160, 181, 250],
         [125, 151, 240],
         [100,  91, 207],
         [ 19,   6, 186],
         [  0,   0, 190],
         [  0,   0, 250],
         [  0,  31, 255],
         [  0,  82, 255],
         [  0, 133, 255],
         [  0, 185, 255],
         [  5, 236, 241],
         [ 47, 255, 200],
         [ 88, 255, 159],
         [129, 255, 117],
         [170, 255,  76],
         [212, 255,  35],
         [254, 237,   0],
         [255, 190,   0],
         [255, 142,   0],
         [255,  95,   0],
         [255,  47,   0],
         [232,   1,   0],
         [174,   0,   0],
         [150,  47,  47],
         [173, 112, 112],
         [227, 205, 205]]

  CLR = np.array(CLR)/255.
  CMP = create_colormap(CLR, nclrs)

  return CMP


def colormap_uv(nclrs=200):
  """
    Colormap for U, V with negative - positive values
    with white in the middle
  """
  from matplotlib import cm
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap

  CLR = [[ 19,   5,  77],
         [ 25,   2, 120],
         [ 39,   5, 179],
         [ 52,   5, 245],
         [ 93,   5, 245],
         [129,   5, 245],
         [177,   5, 245],
         [200,  86, 245],
         [223, 161, 247],
         [241, 217, 250],
         [255, 255, 255],
         [250, 242, 217],
         [250, 217, 110],
         [245, 193,  88],
         [247, 175,  30],
         [235, 138,  42],
         [224,  96,  31],
         [214,  79,  11],
         [168,  59,   5],
         [122,  42,   2],
         [ 92,  32,   2]]

  CLR = np.array(CLR)/255.
  CMP = create_colormap(CLR, nclrs)

  return CMP



def colormap_temp(nclrs=200, clr_ramp=[1,1,1], add_btm=True):
  """
    Colormap for temp
    low S value ramp to clr_ramp
  """
  from matplotlib import cm
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap

  btm = cm.get_cmap('jet',nclrs)
  ixtop  = round(nclrs*0.1)-1
  clrbtm = btm(range(nclrs))
  chbtm  = np.zeros((ixtop,4))
  if add_btm == True:
# Add white at the beginning:
    cxbtm  = clrbtm[0,:]
  else:
# Add white at the top
    ixtop  = round(nclrs*0.1)-1
    ixbtm  = nclrs-ixtop-1
    cxbtm  = clrbtm[ixbtm,:]

  chbtm[:,3] = cxbtm[3]

  for ik in range(3):
    cc0 = clr_ramp[ik]
    chbtm[:,ik]  = np.linspace(cxbtm[ik],cc0,ixtop)

  if add_btm:
    chbtm = np.flip(chbtm, axis=0)
    clrbtm = np.insert(clrbtm,0,chbtm, axis=0)
  else:
    clrbtm[ixbtm+1:nclrs,:] = chbtm

  newclrs = clrbtm
  newcmp  = ListedColormap(newclrs)

  return newcmp

def colormap_ssh(cpos='Oranges',cneg='Blues_r',nclrs=100, clr_ramp=[1,1,1]):
  """
    Create colormaps for showing positive/negative ranges
    clr_ramp - color in the middle of the colormap
    nclrs - # of colors in positve or negative segments
  """
  from matplotlib import cm
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap

  btm = cm.get_cmap(cpos,nclrs)
  ixpos  = round(nclrs*0.1)-1
  clrpos = btm(range(nclrs))
  chpos  = np.zeros((ixpos,4))
  chneg  = np.zeros((ixpos,4))

# Positive colors:
# Add white at the beginning:
  fixclr = clrpos[ixpos,:]
  chpos[:,3] = fixclr[3]

  clr_ramp=[1,1,1]
  for ik in range(3):
    cc0 = clr_ramp[ik]
    chpos[:,ik]  = np.linspace(fixclr[ik],cc0,ixpos)

  chpos = np.flip(chpos, axis=0)
  clrpos[:ixpos,:] = chpos

# Negative ssh range:
# Add white at the top
  clr0    = cm.get_cmap(cneg,nclrs)
  clrneg  = clr0(range(nclrs))
  ixneg   = nclrs - ixpos - 1
  fixnegc = clrneg[ixneg,:]
  chneg[:,3] = fixnegc[3]

  for ik in range(3):
    cc0 = clr_ramp[ik]
    chneg[:,ik]  = np.linspace(fixnegc[ik],cc0,ixpos)

#  chneg = np.flip(chneg, axis=0)
  clrneg[ixneg+1:nclrs,:] = chneg

  newclrs = np.append(clrneg,clrpos, axis=0)
#  newclrs = clrbtm
  newcmp  = ListedColormap(newclrs)

  return newcmp

def colormap_topo1(nclrs=200):
  """
    Colormap for topography
  """
#  import mod_colormaps as mclrs
  CLR = [[59,   43, 105],
         [33,   17, 128],
         [69,   52, 166],
         [39,   10, 204],
         [20,  101, 191],
         [14,  139, 252],
         [45,  171, 201],
         [45,  181, 170],
         [25,  207, 155],
         [89,  194, 120],
         [63,  209,  50],
         [170, 242,  53],
         [242, 247,  89],
         [249, 250, 212]]

  CLR = np.array(CLR)/255.
  CMP = create_colormap(CLR, nclrs)

  return CMP

def colormap_discrete(CLR=[]):
  """
    Create custom colorbar for discrete values
    CLR is 2D list
  """
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap
  if len(CLR) == 0:
    CLR = [[  4,   4,  66],
           [  7,   7, 242],
           [133,   7, 242],
           [242,   7, 223],
           [  7, 195, 242],
           [  7, 242, 105],
           [ 53,  99,  72],
           [240, 224,   7],
           [209, 103,   4],
           [209,  38,   4]]

  CLR   = np.array(CLR)/255.
  nclrs = CLR.shape[0]
  CMP = ListedColormap(CLR)
#  CMP   = create_colormap(CLR, nclrs)

  return CMP

