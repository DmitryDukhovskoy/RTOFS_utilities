# Example to compute RMSD from 2 arrays
# of different sizes X1, X2
# size (X1) > size (X2)
# X1, X2 - observations (or model output) at time instances
# T1, T2 (T1 - more frequent)
# 
# Strategy: subsample X1 to time instances T2
import numpy as np

mu1, sgm1  = 0.5, 2.0
mu2, sgm2  = 0.55, 2.2

dt1 = 10.   # 10-min observations
dt2 = 60.   # 1-hr model output
T1  = np.arange(0.,48.*60.,dt1)
T2  = np.arange(0.,48.*60.,dt2)
N1  = T1.shape[0]
N2  = T2.shape[0]

X1 = np.random.normal(mu1, sgm1, N1) # generate N1 random values, "observations"
X2 = np.random.normal(mu2, sgm2, N2) # N2 random values, "model" 

# Rounding minutes to integers:
iT1 = T1.copy().astype(int)
iT2 = T2.copy().astype(int)

# Searching for indices where time T1 = T2
I = np.array([], dtype=int)
for ii in range(N2):
  i0 = np.where(iT1 == iT2[ii])
  I  = np.append(I,i0)

# Subsample X1 using indices I (i.e. at time = T2):
sX1 = X1[I]

# Check that the arrays are same size:
if sX1.shape[0] == N2:
  print('Subsampled array size = control size = {0}'.\
        format(N2))
else:
# Need to edit code to adjusting the size of X2 to match 
# subsampled array, probably the last value was not subsampled
  raise Exception("Need to adjust size of array X2")

# Now perform RMSD:
rmsd = np.sqrt(1./float(N2) * np.dot((sX1-X2),(sX1-X2)))

print('RMSD = {0}'.format(rmsd))


# Another approach:
sd=(sX1-X2)**2
msd=np.nanmean(sd)
rmsd2=np.sqrt(msd)

print('RMSD 2nd appr= {0}'.format(rmsd))



