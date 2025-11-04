"""
  vertical levels for GFS  MOM6 output
"""
zz = [
   0, 2, 5, 7, 10, 15, 20, 25, 30, 50, 75, 100, 125, 150, 200, 250, 300,
   400, 500, 750, 1000, 1250, 1500, 1750,
   *range(2000, 5001, 500)
]

dzz = np.diff(zz)

zm = dzz*0
zm = -dzz[0]*0.5 + np.cumsum(dzz)


zm = [
  1, 3, 5, 7, 10, 15, 20, 25, 30, 50, 75, 100, 125, 150, 200, 250, 300,
   400, 500, 750, 1000, 1250, 1500,
   *range(2000, 5001, 500)
]  
np.set_printoptions(suppress=True, precision=2)
print(zm)



