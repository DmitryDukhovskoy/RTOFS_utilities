nOB = len(segments)

# Time array for monthly clim:
time_days = np.zeros((12))
for imo in range(1,13):
  dnmb = TMM[imo-1,3]
  time_days[imo-1] = dnmb - dnmb_start + 1

#static = xarray.open_dataset('/work/acr/spear/analysis/ocean_z.static.nc')
grid_spear = xarray.open_dataset('/work/Dmitry.Dukhovskoy/data/SPEAR/ocean_z.static.nc')
icegrid_spear = xarray.open_dataset('/work/Dmitry.Dukhovskoy/data/SPEAR/ice.static.nc')

fconfig = 'config_nep.yaml'
with open(fconfig) as ff:
  config = safe_load(ff)

# Check if mapping indices exist, gmapi:
dirgmapi = config['filesystem']['spear_mom_gmapi']
flgmaph  = f'spear2mom_NEP_OB_gmapi_hpnt.nc'
flgmapu  = f'spear2mom_NEP_OB_gmapi_upnt.nc'
flgmapv  = f'spear2mom_NEP_OB_gmapi_vpnt.nc'
dflgmaph = os.path.join(dirgmapi, flgmaph)
dflgmapu = os.path.join(dirgmapi, flgmapu)
dflgmapv = os.path.join(dirgmapi, flgmapv)
# h-point indices
dsh = xarray.open_dataset(dflgmaph)
# u-point indices
dsu = xarray.open_dataset(dflgmapu)
# v-point indices
dsv = xarray.open_dataset(dflgmapv)

varnm = varplt
nens = len(ENSR)

icc = 0
dsetOB = xarray.Dataset()
for isgm in range(nOB):
  nsgm = isgm+1
  print(f'Processing lon/lat OB segment={nsgm}')
  dset   = mutob.derive_obsegm_lonlat(hgrid, segments, isgm)
  dsetOB = xarray.merge([dsetOB, dset])

for iens in range(nens):
  ens = ENSR[iens]
  spear_dir = config['filesystem']['nep_spear_subset'].\
                   format(year=dv_start[0], ens=ens)

  if varnm == 'thetao' or varnm == 'so':
    # Load monthly SPEAR data subset for NEP
    flnm_spear = f'NEP_spear_{dv_start[0]}{dv_start[1]:02d}.{varnm}.nc'
    ds = mutob.read_spear_output(spear_dir, varnm, flnm_spear, fzint=True)

    for isgm in range(nOB):
      nsgm   = isgm+1
      print(f'Processing {varnm} OB segment={nsgm} ens={ens:02d}')
      INDX   = dsh[f'indx_segm{nsgm:03d}'].data
      JNDX   = dsh[f'jndx_segm{nsgm:03d}'].data
      # Interpolate onto NEP OB supergrid:
      dset   = mutob.derive_obsegm_3D(hgrid, ds, segments, isgm, varnm,
                                      INDX, JNDX, time_steps=time_days)
  #    dset   = xarray.Dataset({f"{varnm}_segment_{nsgm:03d}": darr})
      dsetOB = xarray.merge([dsetOB, dset])
      fldnm_old = f'{varnm}_segment_{nsgm:03d}'
      fldnm_new = f'{varnm}_e{ens:02d}_segment_{nsgm:03d}'
      dsetOB = dsetOB.rename({fldnm_old: fldnm_new})

# UV fields
  elif varnm == 'u' or varnm == 'v':
    # Derive rotation angle, rad and components of the rotation matrix
    # for rotating vectors onto true N/E grid from SPEAR
    r2d = 180./np.pi
    theta_rot, cosrot, sinrot = mutob.get_rotangle(icegrid_spear, fconfig, grid_spear)
    print(f"Rotation angle for SPEAR min/max: {np.min(theta_rot)*r2d:6.2f} " + \
          f"/ {np.max(theta_rot)*r2d:6.2f}")

    # Load monthly SPEAR data subset for NEP
    flnmu_spear = f'NEP_spear_{dv_start[0]}{dv_start[1]:02d}.uo.nc'
    flnmv_spear = f'NEP_spear_{dv_start[0]}{dv_start[1]:02d}.vo.nc'
    ds_uo = mutob.read_spear_output(spear_dir, 'uo', flnmu_spear, fzint=True)
    ds_vo = mutob.read_spear_output(spear_dir, 'vo', flnmv_spear, fzint=True)

    for isgm in range(nOB):
      nsgm   = isgm+1
      print(f'Processing {varnm} OB segment={nsgm}')

      if varnm == 'u':
        INDX  = dsu[f'indx_segm{nsgm:03d}'].data
        JNDX  = dsu[f'jndx_segm{nsgm:03d}'].data
      elif varnm == 'v':
        INDX  = dsv[f'indx_segm{nsgm:03d}'].data
        JNDX  = dsv[f'jndx_segm{nsgm:03d}'].data

      dset   = mutob.derive_obsegm_uv(hgrid, ds_uo, ds_vo, segments, isgm, theta_rot,\
                                      varnm, INDX, JNDX, time_steps=time_days)
      dsetOB = xarray.merge([dsetOB, dset])
      fldnm_old = f'{varnm}_segment_{nsgm:03d}'
      fldnm_new = f'{varnm}_e{ens:02d}_segment_{nsgm:03d}'
      dsetOB = dsetOB.rename({fldnm_old: fldnm_new})

  for segm in [1,2,3,4]:
#    vv  = f"{varnm}_segment_{segm:03d}"
    vdz = f"dz_{varnm}_segment_{segm:03d}"
    dm1 = f'lat_segment_{segm:03d}'
    dm2 = f'lon_segment_{segm:03d}'
    dsetOB[fldnm_new].attrs["coordinates"] = f"{dm1} {dm2}"
    dsetOB[vdz].attrs["coordinates"] = f"{dm1} {dm2}"


f_save = False
if f_save:
  date_init = f'{dv_start[0]}{dv_start[1]:02d}{dv_start[2]:02d}'
  pthoutp = gridfls['MOM6_NEP'][run_name]['pthoutp']
  fobc_out = os.path.join(pthoutp,f'OBCs_spear_mnth_check{date_init}.nc')
  print(f'Saving SPEAR OBCs from {nens} ensemlbe runs ---> {fobc_out}')

  dsetOB.to_netcdf(fobc_out,
                 format='NETCDF3_64BIT',
                 engine='netcdf4',
                 unlimited_dims='time')



