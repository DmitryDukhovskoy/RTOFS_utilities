import matplotlib.pyplot as plt

def plot_prof(Z,T,ctl="",fgn=1):
  plt.ion()
  fig1 = plt.figure(fgn,figsize=(7,8), constrained_layout=False)
  fig1.clf()
  plt.plot(T,Z)
  plt.title(ctl)

  return

def plot_2prof(Z1,T1,Z2,T2,ctl1="",ctl2="",fgn=1):
  plt.ion()
  fig1 = plt.figure(fgn,figsize=(7,8), constrained_layout=False)
  fig1.clf() 
  ax1 = plt.axes([0.1, 0.1, 0.35, 0.8])
  plt.plot(T1,Z1)
  ax1.grid(True)
#  ax1.set_xlabel('cycle/hr')
  ax1.set_title(ctl1)

  ax2 = plt.axes([0.55, 0.1, 0.35, 0.8])
  ax2.plot(T2,Z2)
  ax2.grid(True)
  ax2.set_title(ctl2)

def plot_checkN2Phi(N2,ZZ,Z_phi,Phi,PRB,zbtm,prbnm,krcrd,Rrsb,Cphase,fgn=4):
  """
    Plot N2, N2 smoothed, 1st eigenvector
  """
  import numpy as np
  import mod_utils_rossby as ursb
  import mod_solver as uslv

  N2zF = uslv.runmn(N2,ZZ,mnwnd=3)  # Filter
  N_hrz = np.sqrt(N2)
  N_chr = 3600.*N_hrz   # cycles per hr
  NF_hrz = np.sqrt(N2zF)
  NF_chr = 3600.*NF_hrz

  plt.ion()
  fig4 = plt.figure(fgn,figsize=(8,8), constrained_layout=False)
  fig4.clf()
  ax1 = plt.axes([0.1, 0.1, 0.35, 0.8])
  plt.plot(N_chr,Z_phi)
  plt.plot(NF_chr,Z_phi)
  ax1.grid(True)
  ax1.set_xlabel('cycle/hr')
  uniqnmb = PRB.uniqnmb
  ctl = 'N^2, {1}, U#{2} Zbtm={0:.1f}'.format(zbtm,prbnm,uniqnmb)
  ax1.set_title(ctl)

  ax2 = plt.axes([0.55, 0.1, 0.35, 0.8])
  ax2.plot(Phi,Z_phi)
  ax2.grid(True)
  ctl2 = '$\Phi$(z)'
  ax2.set_title(r'$\Phi$(z)')

  CC  = PRB.CC
  lat = PRB.lat
  lon = PRB.lon
  yr  = PRB.yr
  mo  = PRB.mo
  md  = PRB.md
  gtm = PRB.gtm
  cid = PRB.cid

  Txt = ['Cruise# {1}{0}'.format(cid,CC)]
  Txt.append('lat={0:.2f}N, lon={0:.2f}E'.format(lat,lon))
  Txt.append('{0:02d}/{1:02d}/{2}:{3:.2f}'.format(md,mo,yr,gtm))
  Txt.append('UnqNmb={0}'.format(uniqnmb))
  Txt.append('Seq.Recrd={0}'.format(krcrd))
  Txt.append('Rossby(km)={0:.1f}'.format(Rrsb))
  Txt.append('PhaseSp(m/s)={0:.2f}'.format(Cphase))
  nt=len(Txt)

  x0 = 1
  y0 = 1
  ax3 = plt.axes([0.56,0.12,0.25,0.1])
  ax3.set(xlim=(x0,4.5),ylim=(y0,6))
  for itt in range(nt):
    ssl=Txt[itt]
    x1 = x0
    y1 = y0+itt
    ax3.text(x1,y1,ssl)
  ax3.axis('off')

  return 


