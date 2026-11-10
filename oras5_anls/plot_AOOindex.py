"""
  AOO index from
  https://www2.whoi.edu/site/beaufortgyre/results/arctic-ocean-oscillation-index-aoo/updated-aoo-index/#:~:text=The%20AOO%20is%20non%2Ddimensional,index%20corresponds%20to%20cyclonic%20circulation.&text=1946.

The AOO is non-dimensional: x1.E-6 (ssh difference (cm) between closed SSH isolines in the BG region divided by distance between these isolines (cm). Positive AOO index corresponds to anticyclonic circulation and a negative index corresponds to cyclonic circulation.

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import matplotlib
import matplotlib.colors as colors
PPTHN = '/home/Dmitry.Dukhovskoy/python'
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')

from mod_utils_fig import bottom_text

YRS = []
AOO = []
with open('AOO_index.txt', 'r') as file:
  for line in file:
    # Split each line into components
    parts = line.strip().split()
    if len(parts) == 2:
      year = int(float(parts[0]))
      val = float(parts[1])
      YRS.append(year)
      AOO.append(val)

YRS = np.array(YRS)
AOO = np.array(AOO)
YRSp = YRS+0.5

xtks = [x for x in range(YRS[0],YRS[-1]+1)]
xtk_labels = [str(int(t)) if i % 2 == 0 else '' for i, t in enumerate(xtks)]



plt.ion()
fig1 = plt.figure(1,figsize=(9,8))
plt.clf()
ax1 = plt.axes([0.1, 0.55, 0.85, 0.4])
sttl = f'AOO annual index'

ax1.plot(YRSp, AOO, color=[0.,0.4,0.9])

ax1.grid('on')
ax1.set_xticks(xtks)
ax1.set_xticklabels(xtk_labels)
ax1.set_xlim([1988,2025])
ax1.set_title(sttl)

btx = 'plot_gradBG.py'
bottom_text(btx, pos=[0.08,0.45])




