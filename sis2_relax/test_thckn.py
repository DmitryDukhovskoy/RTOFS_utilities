"""
  Check cell-mean thcikness calculation
  Want m3 (ice) / m2 (cell area)
"""

ai = np.array([0.0286, 0.00938, 0.00933, 0.9412, 0.00067, 0.00012, 0.000011, 0., 0.0005, 0.0023])
if np.sum(ai) > 1.:
  ai = ai / np.sum(ai)

#rho_ice = 905.
#mHice = np.array([69.05, 121.5, 287.5, 985., 1010, 1731, 1854., 2262., 2959., 3277.]) 
#hi = mHice / rho_ice

hi = np.array([0.07629834, 0.13425414, 0.31767956, 1.08839779, 1.1160221 ,\
       1.91270718, 2.04861878, 2.49944751, 3.26961326, 3.62099448])

hi_cell = np.sum(ai*hi)
print(f'hi_cell={hi_cell}')

