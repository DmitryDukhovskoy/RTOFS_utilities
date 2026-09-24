"""
Fit OLS and derive statistical tests 
taking into accoutn a/correlation time scale
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib
import random
from scipy import stats
from scipy.stats import t as tdist

PPTHN = '/home/Dmitry.Dukhovskoy/python'
if len(PPTHN) == 0:
  cwd   = os.getcwd()
  aa    = cwd.split("/")
  nii   = cwd.split("/").index('python')
  PPTHN = '/' + os.path.join(*aa[:nii+1])
sys.path.append(PPTHN + '/MyPython/hycom_utils')
sys.path.append(PPTHN + '/MyPython/draw_map')
sys.path.append(PPTHN + '/MyPython')
sys.path.append(PPTHN + '/MyPython/mom6_utils')

from mod_utils_fig import bottom_text

import numpy as np
import random

def regr_stat(Ydsz, X, B, hac=False):
  # hac - HAC/Newey-West covariance estimate
  # for taking into account a/correlation
  alf0 = B[0]  # intercept
  alf1 = B[1]  # slope

  # Residuals
  resid = Ydsz - X @ B

  n = len(Ydsz)
  p = X.shape[1]

  # OLS residual degrees of freedom
  df_ols = n - p

  # Residual variance
  s2 = np.sum(resid**2) / df_ols

  # Residual ACF
  res_demean = resid - np.mean(resid)

  acf = np.correlate(res_demean, res_demean, mode='full')

  acf = acf[len(res_demean)-1:]
  acf /= acf[0]
  acf_pos = acf[1:]
  i0 = np.where(acf_pos < 0)[0]

  if len(i0) > 0:
    max_lag = i0[0] + 1
  else:
    max_lag = len(acf_pos)

  # Covariance matrix of regression coefficients
  XtX_inv = np.linalg.inv(X.T @ X)

  if hac:
    # Newey-West / HAC covariance
    Z = X * resid[:, None]

    # Lag 0 contribution
    S = Z.T @ Z

    # Add autocovariance contributions
    for lag in range(1, max_lag + 1):
      # Bartlett kernel weight
      weight = 1.0 - lag / (max_lag + 1.0)

      # Lagged covariance
      Gamma = Z[lag:].T @ Z[:-lag]

      S += weight * (Gamma + Gamma.T)

    # HAC covariance matrix
    cov_beta = XtX_inv @ S @ XtX_inv

    print("Statistics using HAC/Newey-West covariance")
    print(f"Maximum HAC lag = {max_lag}")

  else:
    # Conventional OLS covariance
    cov_beta = s2 * XtX_inv
    print("Statistics assuming independent residuals")

  # Standard errors

  se_beta = np.sqrt(np.diag(cov_beta))

  se_alf0 = se_beta[0]
  se_alf1 = se_beta[1]

  # t-statistics
  t_alf0 = alf0 / se_alf0
  t_alf1 = alf1 / se_alf1

  # P-values
  # HAC inference is asymptotic, so use normal distribution.
  # With n ~ 3650 
  if hac:
    p_alf0 = 2.0 * stats.norm.sf(abs(t_alf0))
    p_alf1 = 2.0 * stats.norm.sf(abs(t_alf1))
    zcrit = stats.norm.ppf(0.975)

  else:
    p_alf0 = 2.0 * stats.t.sf(abs(t_alf0), df_ols)
    p_alf1 = 2.0 * stats.t.sf(abs(t_alf1), df_ols)
    zcrit = stats.t.ppf(0.975, df_ols)

  # 95% confidence interval for slope
  alf1_CI = [alf1 - zcrit * se_alf1, alf1 + zcrit * se_alf1]

  # Print statistics
  print(f'alf0     = {alf0:.5f} +/- {se_alf0:.5f}')
  print(f't stat   = {t_alf0:.3f}')
  print(f'p-val    = {p_alf0:.3e}')

  print(f'alf1     = {alf1:.8f} +/- {se_alf1:.8f}')
  print(f't stat   = {t_alf1:.8f}')
  print(f'p-val    = {p_alf1:.3e}')

  print(
      f'95CI     = '
      f'[{alf1_CI[0]:.8f}, {alf1_CI[1]:.8f}]\n'
  )


random.seed(123)
np.random.seed(123)

tau = 365.
A0  = 5.

mu = 0.
sgm_w = 0.6      # white noise std
sgm_a = 0.2      # AR(1) innovation std
phi   = 0.95     # lag-1 correlation

ndays = 365
nyrs  = 10
alfa  = 0.0003

nrecs = ndays*nyrs
t = np.arange(nrecs)

# Seasonal cycle
a_rnd = random.uniform(0.8,1.2)
t_rnd = random.normalvariate(0.,5.)

season = a_rnd*A0*np.sin(2*np.pi*t/(tau+t_rnd))

# AR(1) noise
Ea = np.zeros(nrecs)
for k in range(1,nrecs):
  Ea[k] = phi*Ea[k-1] + np.random.normal(0.,sgm_a)

# White noise
Ew = np.random.normal(0.,sgm_w,nrecs)

# trend
trend = alfa * np.arange(nrecs)

# Full signal:
# seasonal signal + AR1 + white noise + trend
Yts = season + Ea + Ew + trend


# Processing created time 

# Remove seasonality, use detrended time series to get climatology
Ydetr = season + Ea + Ew
Y2d = Ydetr.reshape(nyrs,ndays)
Y_clim = np.mean(Y2d,axis=0)

Y2dsz = Y2d - Y_clim[None,:]
Ydsz  = Y2dsz.reshape(nrecs)

# Add the trend back:
Ydsz += trend

# Check lag-1 correlation:
r_lag1 = np.corrcoef(Ydsz[:-1],Ydsz[1:])[0,1]
print(f"Estimate lag1 acorr: {r_lag1:.3f}, prescribed acorr for AR1: {phi:.3f}")

# Estimate trend:

X = np.column_stack((np.ones(len(t)), t))

# OLS coefficients
B = np.linalg.lstsq(X, Ydsz, rcond=None)[0]
alf0 = B[0]
alf1 = B[1]

print(f"Estimated slope: {alf1}, true slope: {alfa}")

# Statistics ignoring a/corrleation:
regr_stat(Ydsz, X, B)

# Statistics accounting for a/correlation
regr_stat(Ydsz, X, B, hac=True)



