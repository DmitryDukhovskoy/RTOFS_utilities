"""
  Colormaps
"""
import numpy as np
from copy import copy
import matplotlib as mtplt
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

def minmax_clrmap(dmm, pmin=10, pmax=90, cpnt=0.01, fsym=False):
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

def create_colormap(CLR, Ncmp, cmp_obj=True, add_alpha=False):
  """
    Mix main colors in CLR by linear interpolation to create a smoother colormap.
    CLR      : array of base RGB colors, shape (n,3) or (n,4), 4th position - alfa, transparency
    Ncmp     : desired number of colors in final colormap
    cmp_obj  : True --> return ListedColormap, False --> return RGB array
    add_alpha: True --> add A to RGB if missing, False --> keep RGB
  """
  import numpy as np
  from matplotlib.colors import ListedColormap

  # Convert to array
  CLR = np.array(CLR)
  nClr = CLR.shape[0]
  nClr, nCh = CLR.shape  # nCh = 3 (RGB) or 4 (RGBA)

  # Ensure RGBA (add alpha if missing)
  if nCh == 3 and add_alpha:
    CLR = np.hstack([CLR, np.ones((nClr, 1))])
    nCh = CLR.shape[1]

  # If fewer requested colors than base colors --> return as-is
  if Ncmp <= nClr:
    print('create_colormap:')
    print(f'Specified N of colors {Ncmp} <= N of Main Colors {nClr}')
    print(' Colormap not changed')
    return ListedColormap(CLR) if cmp_obj else CLR

  # Smooth interpolation across the entire color sequence
  #
  # Positions of base colors (0..1)
  base_pos = np.linspace(0, 1, nClr)
  # Desired positions for output colors
  new_pos = np.linspace(0, 1, Ncmp)

  # Allocate full colormap
  newCLR = np.zeros((Ncmp, nCh))

  # Interpolate each channel (R,G,B,A)
  for k in range(nCh):
    newCLR[:,k] = np.interp(new_pos, base_pos, CLR[:,k])

  # Return either the colormap object or raw array
  if cmp_obj:
    return ListedColormap(newCLR)
  else:
    return newCLR


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

def colormap_ice_thkn():
  """
   Prepare colormap for sea ice thickness
   Creates ListedColormap object
   to get colors:
   cmpice.colors
   cmpice.N - # of colors
  
  """
#  import mod_colormaps as mclrs
  CLR = [[255, 255, 255],
         [204, 204, 255],
         [153, 51,  255],
         [102, 102, 255],
         [51,  51,  255],
         [0,   0,   204],
         [0,   204, 102],
         [0,   255,   0],
         [153, 255,  51],
         [255, 255,   0],
         [255, 255, 204],
         [255, 204, 153],
         [255, 128,   0],
         [255, 102, 102],
         [150, 0,     0]]
  
  CLR = np.array(CLR)/255.
  CLR = smooth_colors(CLR, smooth_wnd=0.1, nsmooth=1)
#  CLR = np.flip(CLR, axis=0)
  CLR[0,:]=[1,1,1]  # keep white
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

def colormap_haline2(nclrs=200, end_clr=None, start_clr=[0, 0, 153]):
  """
    Colormap for salinity
    optional: specify start / end colors
     end_clr = [R,G,B]
     start_clr = [R,G,B]
  """
#  import mod_colormaps as mclrs
  CLR = [[ 0,    0, 153],
         [69,   52, 166],
         [124, 103, 226],
         [102, 125, 170],
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

  CLR = np.array(CLR, dtype=float)/255.
  if start_clr is not None and len(start_clr) > 0:
    start_clr = np.array(start_clr, dtype=float)
    if np.max(start_clr) > 1.:
      start_clr /= 255.

    CLR[0, :] = start_clr

  if end_clr is not None and len(end_clr) > 0:
    end_clr = np.array(end_clr, dtype=float)
    if np.max(end_clr) > 1.:
      end_clr /= 255.

    CLR[-1, :] = end_clr

  CMP = create_colormap(CLR, nclrs)

  return CMP

def colormap_salin2(nclrs=200):
  """
    Colormap for salinity
    bright color palett
  """
#  import mod_colormaps as mclrs
  CLR = [[255, 255, 255],
         [255, 255, 255],
         [255, 200, 255],
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
         [255, 120, 50],
         [255, 50, 0],
         [100,  0,  0]]

  CLR = np.array(CLR)/255.
  CLR = smooth_colors(CLR, smooth_wnd=0.2, nsmooth=1)
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
         [185, 150, 205],
         [130, 120, 215],
         [ 80,  90, 200],
         [ 60, 110, 220],
         [ 40, 140, 235],
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

def addendclr_colormap(clrmp_name, clr_ramp, nclrs=200, ramp_start=True, \
                       nramp=0.1, cmp_obj=True):
  """
    Modify existing colormap by adding color=clr_ramp 
    at the start (ramp_start) or end of the colormap
    nramp = ramp rate, bigger nramp - slower ramping
  """
#  from matplotlib.colors import ListedColormap, LinearSegmentedColormap
  cmap = mtplt.colormaps.get_cmap(clrmp_name)
  X    = np.linspace(0,1,nclrs)
  CLR  = cmap(X)[:,:3]
  ixtop = round(nclrs*nramp)-1
  if ixtop < 1: ixtop = 1
  chramp = np.zeros((ixtop, 3))  
  if ramp_start:
    clr_end = CLR[ixtop,:]   
    for ik in range(3):
      cc0 = clr_ramp[ik]
      chramp[:,ik] = np.linspace(cc0, clr_end[ik], ixtop)     
    CLR[:ixtop,:] = chramp
  else:
    clr_end = CLR[-ixtop,:]
    for ik in range(3):
      cc0 = clr_ramp[ik]
      chramp[:,ik] = np.linspace(clr_end[ik], cc0, ixtop)
    CLR[-ixtop:,:] = chramp

  if cmp_obj:
    CMP = ListedColormap(CLR)
  else:
    CMP = np.array(CLR)

  return CMP

def smooth_colors(CLR, smoothS=0., nsmooth=1, smooth_wnd=0.1):
  """
    CLR - 3D color indices
    Smooth colors, smoothing window = 0.1*smoothing range
    Smoothing range =  [smoothS*nclrs : end]
    # of smoothing = nsmooth (box filtering) 
  """
  nclrs  = CLR.shape[0]
  iS     = int(np.floor(smoothS*nclrs))
  iE     = nclrs
  sm_int = iE-iS
  dii   = int(smooth_wnd*sm_int/2)  # half smooth. window
  if dii<1: dii=1
  for jff in range(nsmooth):
    CLR_sm = CLR.copy()
    for kk in range(iS, iE):
      ix1 = kk-dii
      if ix1<0: ix1=0
      ix2 = kk+dii
      if ix2>nclrs-1: ix2=nclrs-1 
      B = CLR[ix1:ix2+1,:]
      CLR_sm[kk,:] = np.mean(B, axis=0)

    CLR = CLR_sm.copy()

  return CLR

def colormap_temp(nclrs=200, clr_ramp=[1,1,1], add_btm=True, skip_dark=0.12):
  """
    Colormap for temp
    low S value ramp to clr_ramp

    nclrs : Total number of colors in the colormap.
    clr_ramp :  RGB color at the bottom of the ramp (values 0–1). Default is white [1,1,1].
    add_btm :   If True, ramp is added at the beginning (low end); if False, at the top (high end).
    skip_dark_frac :   Fraction of the original 'jet' colormap to skip at the low end to remove very dark blues.
        Typical value: 0.1–0.15.

  """
  from matplotlib import cm
  from matplotlib.colors import ListedColormap 


  jet = cm.get_cmap('jet')
  # Skip the darkest part of jet to avoid white -> dark blue jump
  start_frac = skip_dark
  clrbtm = jet(np.linspace(start_frac, 1.0, nclrs))

  # Size of the ramp section
  ramp_size = int(round(nclrs * 0.1))  # 10% of colormap by default

  # Ramp colors:
  ramp = np.zeros((ramp_size, 4))
  if add_btm:
    # Add ramp at the beginning
    start_color = clrbtm[0, :3]  # first color after truncation
  else:
    # Add ramp at the top
    start_color = clrbtm[-1, :3]  # last color of jet

  # Interpolate RGB from start_color to base_color
  for i in range(3):
    ramp[:, i] = np.linspace(start_color[i], clr_ramp[i], ramp_size)
  ramp[:, 3] = 1.0  # alpha channel

  # Insert ramp into colormap
  if add_btm:
    ramp = np.flipud(ramp)
    clrbtm = np.vstack((ramp, clrbtm))
  else:
    clrbtm = np.vstack((clrbtm, ramp))

  newcmp = ListedColormap(clrbtm)
  return newcmp

def colormap_ssh(cpos='Oranges',cneg='Blues_r',nclrs=100, clr_ramp=[1,1,1]):
  """
    Create colormaps for showing positive/negative ranges
    clr_ramp - color in the middle of the colormap
    nclrs - # of colors in positve or negative segments
  """
  from matplotlib import cm

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

def colormap_discrete(CLR=[], cmp_obj=True):
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
#  CMP = ListedColormap(CLR)
#  CMP   = create_colormap(CLR, nclrs)
  if cmp_obj:
    CMP = ListedColormap(CLR)
  else:
    CMP = np.array(CLR)

  return CMP

def clrmp_lmask(nclrs=2,clr_land=[0.3,0.3,0.3]):
  """
    Create colormap for land mask with 2 colors
  """

  from matplotlib import cm
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap

  r, g, b = clr_land[0:3]
  clrs   = cm.get_cmap('GnBu_r',nclrs)
  newclr = clrs(range(nclrs))
  newclr[0,:] = [1, 1, 1, 1]
  newclr[1,:] = [r, g, b, 1]

  newcmp  = ListedColormap(newclr)

  return newcmp

def colormap_salin(nclrs=200, clr_ramp=[1,1,1]):
  """
    Colormap for salinity
    low S value ramp to clr_ramp
  """
  from matplotlib import cm
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap
  import mod_colormaps as mclrs

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
  CMP = mclrs.create_colormap(CLR, ixtop, cmp_obj=False)
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

def positive_negative(nclrs=200, cname='tgv', cmp_obj=True, neutr=None):
  """
    Several diverging colormaps 
    for showing positive-negative values
    neutr    : specifies neutral color other than default
    
    Colormaps:
    tgv      :  teal - grey - violet
    gwm      : green - white - purple (colorblind safe), similar to tgv
    bwb      : brown - white - blue
    gyp      : green - yellow - purple
    cbo      : cayn - black - orange (maximum contrast)
    bwo      : navy - white - orange
    tbbo     : Teal–Blue - Light - Burnt Orange

    cmp_obj  : True --> return ListedColormap, False --> return RGB array    
  """
  match cname:
    case 'tgv':
      clr_btm = [[  0, 100,  95],
                 [ 80, 160, 150]]
      clr_mid =  [230, 230, 230]
      clr_top = [[150, 110, 170],
                 [ 90,  60, 120]]
    case 'gwm':
      clr_btm = [[  0, 120,  60],
                 [120, 200, 140]]
      clr_mid =  [255, 255, 255]
      clr_top = [[200, 120, 200],
                 [130,  40, 130]]
    case 'bwb':
      clr_btm = [[120,  70,  40],
                 [200, 150, 110]]
      clr_mid = [245, 245, 245]
      clr_top = [[120, 160, 210],
                   [ 40,  80, 160]]
    case 'gyp':
      clr_btm = [[  0, 100,  40],
                 [120, 180,  90]]
      clr_mid = [255, 245, 200]
      clr_top = [[170, 120, 200],
                 [100,  50, 150]]
    case 'cbo':
      clr_btm = [[230, 150,  60],
                 [180,  80,   0]]
      clr_mid = [ 30,  30,  30]
      clr_top = [[  0, 180, 180],
                 [100, 230, 230]]
    case 'bwo':
      clr_btm = [[200,  90,  30],
                 [240, 180, 120]]
      clr_mid = [245, 245, 245]
      clr_top = [[ 80, 120, 190],
                 [ 20,  40, 120]]
    case 'tbbo':
      clr_btm = [[210, 110,  40],
                 [250, 190, 120]]
      clr_mid = [250, 250, 250]
      clr_top = [[ 80, 150, 210],
                 [  0,  90, 150]]
    case _:
      raise ValueError(f"unrecognized colormap option {cname}")

  
  if neutr is not None:
    clr_mid = neutr

  n_half = nclrs // 2
  # Interpolate bottom and top halves
  clr1 = create_colormap(clr_btm + [clr_mid], n_half + 1, cmp_obj=False)
  clr2 = create_colormap([clr_mid] + clr_top, n_half + 1, cmp_obj=False)

  # Remove duplicated middle color
  CLR = np.vstack([clr1[:-1], clr_mid, clr2[1:]])
  CLR = np.flipud(CLR)

  # Normalize to [0,1]
  CLR = CLR / 255.0

  return ListedColormap(CLR) if cmp_obj else CLR


def colormap_temp_old(nclrs=200, clr_ramp=[1,1,1], add_btm=True):
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

def colormap_posneg_uneven(CLRS, nclrs=200):
  """
    Create colormaps for showing positive/negative ranges
    with uneqeual numbers of positive and negative colors, 

   CLRS - specify main colors including the 0-color (white, e.g.)
   if not, use default
    nclrs = desired total # of colors
  """
  from matplotlib import cm
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap

  if len(CLRS) == 0:
    CLRS = [[0.6, 0.02, 0.6],
            [0.2, 0.38, 1],
            [1, 1, 1],
            [0., 0.8, 0.8],
            [0.4, 0.8, 0],
            [1, 1, 0.5],
            [1, 0.8, 0.6],
            [1, 0.6, 0],
            [0.7, 0.1, 0.1]]

  clrmp = create_colormap(CLRS, nclrs)

  return clrmp

def minmax_clrmap(dmm, pmin=10, pmax=90, cpnt=0.01, fsym=False):
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

def colormap_cold(CLRMP=['Purples_r','BuPu','Blues_r'], clrE=[1,1,1], nclrs=50):
  """
    Colormap for negative  values
    Specify colormaps to combine
    start with coldest
    nclrs = # of colors in each colormap
    clrE = color to ramp the end of the new colormap
    colormaps are smoothed at the edges
  """
  nmp = len(CLRMP)
  sclr1 = [0.1,0,0.2]  # start color of coldest range
  eclr1 = [0., 0., 0.1]  # end color for coldest range
  for ii in range(nmp):
    clrmp_name = CLRMP[ii]
    if ii == 0:
#      CLR = addendclr_colormap(clrmp_name, clrS, nclrs=nclrs, cmp_obj=False)
      cmap = mtplt.colormaps.get_cmap(clrmp_name)
      X    = np.linspace(0,1,nclrs)
      CLR  = cmap(X)[:,:3]
      for ipp in range(1,7):
        CLR[ipp-1,:] = sclr1
        CLR[-ipp,:] = eclr1  # modify end color
      CLR = smooth_colors(CLR, smooth_wnd=0.2, nsmooth=1) 
    else:
      clr_next = addendclr_colormap(clrmp_name, CLR[-1,:], \
                   nclrs=nclrs, nramp=0.3, cmp_obj=False)   
      CLR = np.append(CLR,clr_next, axis=0)
 
  for ipp in range(1,5):
    CLR[-ipp,:] = clrE
    CLR[ipp-1,:] = sclr1
#  iS  = int(len(CLR) - 0.4*nclrs)
  iS = 0
  CLR = smooth_colors(CLR, smoothS=iS, smooth_wnd=0.15)
#  CLR[-1,:] = clrE
#  CLR = smooth_colors(CLR, smoothS=iS, smooth_wnd=0.25)
  CMP = ListedColormap(CLR) 
  return CMP

def colormap_warm(CLRMP=['summer','Wistia','gist_heat_r'], clrS=[1,1,1], nclrs=50):
  """
    Colormap for positive  values
    Specify colormaps to combine
    nclrs = # of colors in each colormap
    clrE = color to ramp the end of the new colormap
    colormaps are smoothed at the edges
  """
  nmp = len(CLRMP)
  sclr1 = [1.,1.,1.]  # start color of warm
  clrE  = [0.5, 0.2, 0] # warmest color
#  eclr1 = [0., 0., 0.1]  # end color for warm range

  for ii in range(nmp):
    clrmp_name = CLRMP[ii]
#    cmap = mtplt.colormaps.get_cmap(clrmp_name)
#    X    = np.linspace(0,1,nclrs)
    if ii == 0:
      cmap = mtplt.colormaps.get_cmap(clrmp_name)
      X    = np.linspace(0,1,nclrs)
      CLR  = cmap(X)[:,:3]
      for ipp in range(7):
        CLR[ipp,:] = sclr1
#        CLR[-ipp,:] = eclr1  # modify end color
      CLR = smooth_colors(CLR, smooth_wnd=0.2, nsmooth=1)
    else:
      clr_next = addendclr_colormap(clrmp_name, CLR[-1,:], \
                   nclrs=nclrs, nramp=0.3, cmp_obj=False)   
      CLR = np.append(CLR,clr_next, axis=0)
   
#  for ipp in range(1,5):
#    CLR[-ipp,:] = clrE
  CLR[-1,:] = clrE

  for ipp in range(7):
    CLR[ipp,:] = sclr1
#  iS  = int(len(CLR) - 0.4*nclrs)
  iS = 0
  CLR = smooth_colors(CLR, smoothS=iS, smooth_wnd=0.15)
  CLR[0,:] = clrS
 
  CMP = ListedColormap(CLR) 
  return CMP

def colormap_albedo(Ncmp=200):
  """
    The progression goes from 
     dark gray/black --> blue-gray --> ice blue --> white
     which visually matches the concept of increasing albedo
  """
  CLR_albedo = np.array([
      [0.05, 0.05, 0.05],
      [0.15, 0.20, 0.25],
      [0.25, 0.35, 0.45],
      [0.35, 0.50, 0.60],
      [0.50, 0.65, 0.80],
      [0.65, 0.78, 0.90],
      [0.75, 0.85, 0.95],
      [0.85, 0.90, 0.98],
      [0.93, 0.96, 1.00],
      [1.00, 1.00, 1.00],
  ])

  CMP = create_colormap(CLR_albedo, Ncmp, cmp_obj=True, add_alpha=False)

  return CMP

def colormap_blue_yellow(Ncmp=200):
  """
    perceptually uniform and goes from cold (blue) to yellow
  """

  CLR_thermal = np.array([
    [0.0,   0.043, 0.208],
    [0.0,   0.115, 0.303],
    [0.0,   0.184, 0.395],
    [0.016, 0.253, 0.482],
    [0.047, 0.319, 0.565],
    [0.090, 0.385, 0.645],
    [0.145, 0.447, 0.720],
    [0.211, 0.509, 0.791],
    [0.288, 0.567, 0.857],
    [0.374, 0.622, 0.917],
    [0.469, 0.673, 0.969],
    [0.571, 0.719, 1.000],
    [0.672, 0.762, 0.973],
    [0.765, 0.800, 0.934],
    [0.848, 0.835, 0.889],
    [0.918, 0.867, 0.839],
    [0.969, 0.895, 0.785],
    [0.992, 0.922, 0.723],
    [0.996, 0.941, 0.654],
    [1.0,   0.960, 0.580]
  ])

  CMP = create_colormap(CLR_thermal, Ncmp, cmp_obj=True, add_alpha=False)

  return CMP
 
def colormap_cold_warm(Ncmp=200, ins_white=False, n_white=3):
  """
    Colormap for showing temperatures etc.
    ins_white = True - white in the middle for showing -/+ values
    n_white - width of the color region
  """
  CLR_thermal = np.array([
  [0.2, 0.0, 0.5],    # dark purple
  [0.3, 0.0, 0.6],    # purple
  [0.1, 0.2, 0.7],    # indigo
  [0.0, 0.3, 0.7],    # blue
  [0.0, 0.4, 0.8],    # blue-cyan
  [0.0, 0.5, 0.85],   # cyan-blue
  [0.0, 0.6, 0.75],   # cyan-teal
  [0.0, 0.7, 0.6],    # teal-green
  [0.2, 0.8, 0.5],    # green-teal
  [0.4, 0.85, 0.3],   # green
  [0.55, 0.9, 0.1],   # yellow-green
  [0.7, 0.85, 0.0],   # yellow
  [0.85, 0.7, 0.0],   # orange-yellow
  [0.95, 0.5, 0.0],   # orange
  [0.98, 0.3, 0.0],   # red-orange
  [1.0, 0.15, 0.0],   # red
  [0.9, 0.0, 0.0],    # dark red
  [0.7, 0.0, 0.0]     # very dark red
  ])

  nclrs, nCh = CLR_thermal.shape
  if ins_white:
    ins0 = nclrs // 2
    white = np.ones((n_white, 3))
    CLR_thermal = np.insert(
        CLR_thermal,
        ins0,
        white,
        axis=0
    )

  CMP = create_colormap(CLR_thermal, Ncmp, cmp_obj=True, add_alpha=False)

  return CMP

def colormap_difference_negpos(Ncmp=200, n_white=4):
  """
    Colormap for showing differences with positive and negative values centered around 0
    or (if no white) to show any other field from dark purple -> blue -> yellow -> dark red

    For negative positive:
    transition - white
    n_white - width of the color region

    n_white = 0 - no white
  """

  CLR_thermal = np.array([
      [0.20, 0.00, 0.50],   # dark purple
      [0.30, 0.00, 0.60],   # purple
      [0.10, 0.20, 0.70],   # indigo
      [0.00, 0.30, 0.70],   # blue
      [0.00, 0.45, 0.85],   # blue-cyan
      [0.00, 0.60, 0.85],   # cyan
      [0.20, 0.75, 0.85],   # light cyan
      [0.65, 0.85, 0.90],   # pale blue
      [1.00, 0.95, 0.60],   # pale yellow
      [0.95, 0.85, 0.10],   # yellow
      [0.90, 0.65, 0.00],   # yellow-orange
      [0.95, 0.45, 0.00],   # orange
      [0.98, 0.25, 0.00],   # red-orange
      [1.00, 0.10, 0.00],   # red
      [0.80, 0.00, 0.00],   # dark red
  ])
   
  # Explicitly define the negative/positive transition
  neg = CLR_thermal[:8]
  pos = CLR_thermal[8:]

  # Number of bins on each side
  n_side = (Ncmp - n_white) // 2

  # Interpolate each side independently
  neg_cmap = create_colormap(neg, n_side, cmp_obj=True, add_alpha=False)
  pos_cmap = create_colormap(
      pos, Ncmp - n_white - n_side,
      cmp_obj=True, add_alpha=False
  )

  neg_colors = neg_cmap(np.linspace(0, 1, n_side))[:, :3]
  pos_colors = pos_cmap(np.linspace(0, 1, Ncmp - n_white - n_side))[:, :3]

  if n_white > 0:
    white = np.ones((n_white, 3))
    colors = np.vstack([
        neg_colors,
        white,
        pos_colors
    ])
  else:
    colors = np.vstack([
        neg_colors,
        pos_colors
    ])      


  CMP = create_colormap(
      colors,
      Ncmp,
      cmp_obj=True,
      add_alpha=False
  )

  return CMP

def colormap_temperature_coldwarm(Ncmp=200):
  """
    Colormap for showing scalar fields from low to high
    with contrast colormap: dark purple --> blue --> green --> orange --> red
  """
  CLR_thermal = np.array([
      [0.12, 0.00, 0.35],   # very dark purple
      [0.20, 0.00, 0.50],   # dark purple
      [0.32, 0.00, 0.60],   # purple
      [0.25, 0.08, 0.68],   # violet
      [0.10, 0.20, 0.70],   # indigo
      [0.00, 0.30, 0.70],   # blue
      [0.00, 0.42, 0.78],   # blue
      [0.00, 0.55, 0.82],   # blue-cyan
      [0.00, 0.65, 0.80],   # cyan
      [0.10, 0.72, 0.75],   # turquoise
      [0.30, 0.78, 0.70],   # teal-green
      [0.50, 0.83, 0.70],   # green
      [0.68, 0.86, 0.72],   # pale green
      [0.82, 0.89, 0.70],   # yellow-green
      [0.95, 0.94, 0.65],   # pale yellow
      [1.00, 0.90, 0.40],   # yellow
      [0.95, 0.78, 0.15],   # yellow-orange
      [0.90, 0.62, 0.00],   # orange
      [0.95, 0.42, 0.00],   # orange-red
      [0.98, 0.22, 0.00],   # red-orange
      [1.00, 0.08, 0.00],   # red
      [0.80, 0.00, 0.00],   # dark red
  ])

  CMP = create_colormap( 
      CLR_thermal,
      Ncmp,
      cmp_obj=True, 
      add_alpha=False
  )      
         
  return CMP




