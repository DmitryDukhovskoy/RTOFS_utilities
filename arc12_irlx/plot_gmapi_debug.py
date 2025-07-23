"""
  Quick plot of PIOMAS grid and problematic grid point
  to debug find_PIOMAS_to_ARC12_gmapi.py

  For interactive debugging, LON, LAT< etc are already in python memory

  For getting closest pnts and box vertices - run interactively mod_regmom.py: find_gridpnts_box
"""
import mod_misc1 as mmisc1
import mod_bilinear as mblnr

dhstep=1

plt.ion()


fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
ax1.plot(LON,LAT,'.',color=[0.9,0.9,0.9])
ax1.plot(x0,y0,'r*')

# Closest pnt: xv1,yv1 - geogr coord, iv1,jv1 - indices in PIOMAS grid
ax1.plot(xv1,yv1,'ro')
plot_p1_m1(iv1,jv1,LON,LAT,ax1)


# Plot in index space:
jpdim, ipdim = LON.shape
IP, JP = np.meshgrid(np.arange(ipdim), np.arange(jpdim))
ax1.cla()
ax1.plot(IP,JP,'.',color=[0.9,0.9,0.9])
ax1.plot(iv1,jv1,'ro')
ax1.contour(IP, JP, LON, levels=[xv1], linestyles='solid', colors=[(1, 0, 0)])
ax1.contour(IP, JP, LAT, levels=[yv1], linestyles='solid', colors=[(1, 0, 0)])

plot_p1_m1_ispace(iv1,jv1,IP,JP,ax1)


# Try a simple approach first:
INp = False
if not INp:
  xv1 = LON[jv1,iv1]
  yv1 = LAT[jv1,iv1]
  IV,JV,INp = find_box_include([x0,y0],[iv1,jv1],LON,LAT, plot_bxs=True, ax1=ax1)


def plot_p1_m1(iv1,jv1,LON,LAT,ax1):
  xp = LON[jv1,iv1-1]
  yp = LAT[jv1,iv1-1]
  ax1.plot(xp,yp,'o')
  ax1.text(xp,yp,'i-1')

  xp = LON[jv1,iv1+1]
  yp = LAT[jv1,iv1+1]
  ax1.plot(xp,yp,'o')
  ax1.text(xp,yp,'i+1')

  xp = LON[jv1-1,iv1]
  yp = LAT[jv1-1,iv1]
  ax1.plot(xp,yp,'o')
  ax1.text(xp,yp,'j-1')

  xp = LON[jv1+1,iv1]
  yp = LAT[jv1+1,iv1]
  ax1.plot(xp,yp,'o')
  ax1.text(xp,yp,'j+1')

  ax1.plot([LON[jv1,iv1-1],LON[jv1,iv1+1]],[LAT[jv1,iv1-1],LAT[jv1,iv1+1]],'-')
  ax1.plot([LON[jv1-1,iv1],LON[jv1+1,iv1]],[LAT[jv1-1,iv1],LAT[jv1+1,iv1]],'-')

def plot_p1_m1_ispace(iv1,jv1,IP,JP,ax1):
  ax1.plot(iv1-1,jv1,'o')
  ax1.text(iv1-1,jv1,'i-1')

  ax1.plot(iv1+1,jv1,'o')
  ax1.text(iv1+1,jv1,'i+1')

  ax1.plot(iv1,jv1+1,'o')
  ax1.text(iv1,jv1+1,'j+1')

  ax1.plot(iv1,jv1-1,'o')
  ax1.text(iv1,jv1-1,'j-1')

  ax1.plot([iv1-1,iv1+1],[jv1,jv1],'-')
  ax1.plot([iv1,iv1],[jv1-1,jv1+1],'-')


def find_box_include(XY0,IJ1,LON,LAT, plot_bxs=False, ax1=[]):
  """
    Find a grid cell that encloses a pnt XY0
    given the first nearst vertex 
  """
  import mod_misc1 as mmisc1
  import mod_bilinear as mblnr

  x0, y0   = XY0
  iv1, jv1 = IJ1
  BX = np.array([
      [[iv1,   jv1],
       [iv1,   jv1-1],
       [iv1-1, jv1-1],
       [iv1-1, jv1]],

      [[iv1,   jv1],
       [iv1+1, jv1],
       [iv1+1, jv1-1],
       [iv1,   jv1-1]],

      [[iv1,   jv1],
       [iv1,   jv1+1],
       [iv1+1, jv1+1],
       [iv1+1, jv1]],

      [[iv1,   jv1],
       [iv1-1, jv1],
       [iv1-1, jv1+1],
       [iv1,   jv1+1]]
  ])

  for ibox in range(4):
    IV = BX[ibox,:,0]
    JV = BX[ibox,:,1]
    XX = LON[JV,IV]
    YY = LAT[JV,IV]
    XXc = 0.25*np.sum(XX)
    YYc = 0.25*np.sum(YY)
    XV, YV  = mblnr.lonlat2xy_wrtX0(XX, YY, XXc, YYc)
    x0c, y0c = mblnr.lonlat2xy_pnt(x0,y0, XXc, YYc)
    INp     = mmisc1.inpolygon_1pnt(x0c, y0c, XV, YV)

    if plot_bxs:
      ax1.plot(XX,YY,'.-')
      ax1.plot(XXc,YYc,'x')
      ax1.text(XXc,YYc,f'Box {ibox}')

    if INp: break

  if not INp:
    IV = []
    JV = []

  return IV,JV,INp






img = ax1.pcolormesh(xR, yR, RLXHR, cmap=clrmp, vmin=rmin, vmax=rmax)
#  img = ax1.pcolormesh(RLXHR, cmap=clrmp)

ax1.set_title('Relaxation time, hrs')

ax2 = fig1.add_axes([ax1.get_position().x1+0.025, ax1.get_position().y0,
                   0.02, ax1.get_position().height])
# extend: min, max, both
clb = plt.colorbar(img, cax=ax2, orientation='vertical', extend='both')
ax2.yaxis.set_ticks(list(np.linspace(rmin,rmax,11)))
ax2.set_yticklabels(ax2.get_yticks())
ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
clb.ax.set_yticklabels(["{:.1f}".format(i) for i in clb.get_ticks()], fontsize=10)
clb.ax.tick_params(direction='in', length=12)

ax1.set_title('Searching for a box enclosing MOM6 grid point, PIOMAS grid')
btx = 'plot_gmapi_debug.py'
bottom_text(btx, pos=[0.2, 0.01])


