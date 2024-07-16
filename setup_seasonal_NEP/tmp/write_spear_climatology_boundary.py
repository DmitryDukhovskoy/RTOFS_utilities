import datetime as dt
import numpy as np
from pathlib import Path
import xarray
import os

import sys
sys.path.append('./seasonal-workflow')
from boundary import Segment
from yaml import safe_load

indir = Path('/work/acr/spear/processed/ensmean')
outdir = Path('/work/acr/spear/climatology')

fyaml = 'pypaths_gfdlpub.yaml'
with open(fyaml) as ff:
  gridfls = safe_load(ff)

pthtopo = gridfls['MOM6_NEP']['test']['pthgrid']
fgrid   = gridfls['MOM6_NEP']['test']['fgrid']

hgrid = xarray.open_dataset(os.path.join(pthtopo,fgrid))
segments = [
    Segment(1, 'north', hgrid, output_dir=outdir),
    Segment(2, 'east',  hgrid, output_dir=outdir),
    Segment(3, 'south', hgrid, output_dir=outdir),
    Segment(4, 'west',  hgrid, output_dir=outdir)
]

static = xarray.open_dataset('/work/acr/spear/analysis/ocean_z.static.nc')


def load_var(var):
    if var == 'ssh':
        spear_var = xarray.open_dataset(indir / f'ocean.monthly_mean.ensmean.{var}.nc')
    else:
        spear_var = (
            xarray.open_dataset(indir / f'ocean_z.monthly_mean.ensmean.{var}.nc')
            .rename({'z_l': 'z'})
        )
    # Make sure missing/fill value was properly used
    assert np.abs(spear_var).max() < 1e5
    return spear_var[var]


def compute_climatology(da):
    climo = (
        da
        .sel(start=slice('1993', '2019')) # note hard coded climatology period
        .groupby('start.month')
        .mean('start')
    )
    return climo


def add_valid_time(ds):
    months = ds.lead.values.astype(int) + mon
    dates = [dt.datetime(1993, m, 1) if m <= 12 else dt.datetime(1994, m-12, 1) for m in months]
    ds['time'] = (('lead', ), dates)
    ds = ds.swap_dims({'lead': 'time'})
    return ds


def add_coords(da):
    xcoord = [str(d) for d in da.coords if d in ['xh', 'xq']]
    ycoord = [str(d) for d in da.coords if d in ['yh', 'yq']]
    assert len(xcoord) == len(ycoord) == 1
    xc = xcoord[0]
    yc = ycoord[0]
    coord_tuple = (yc, xc)
    coord_sel = {xc: da[xc], yc: da[yc]}
    lat = next(static[v] for v in ['geolat', 'geolat_u', 'geolat_v'] if xc in static[v].coords and yc in static[v].coords)
    lon = next(static[v] for v in ['geolon', 'geolon_u', 'geolon_v'] if xc in static[v].coords and yc in static[v].coords)
    da['lat'] = (coord_tuple, lat.sel(**coord_sel).data)
    da['lon'] = (coord_tuple, lon.sel(**coord_sel).data)
    return da


if __name__ == '__main__':
    for var in ['ssh', 'thetao', 'so']:
        print(var)
        ds = load_var(var)
        climo = compute_climatology(ds)
        climo = add_coords(climo)
        # this wont work if using daily SSH
        for mon in [3, 6, 9, 12]:
            print(f'  month {mon:02d}')
            mon_data = climo.sel(month=mon)
            mon_data = add_valid_time(mon_data)
            for seg in segments:
                print('    ' + seg.segstr)
                regridded = seg.regrid_tracer(
                    mon_data,
                    regrid_suffix='spear_tracer',
                    write=False
                )
                seg.to_netcdf(regridded, f'spear_{var}_i{mon:02d}_climo')
    print('uv')
    u_climo = add_coords(compute_climatology(load_var('uo')))
    v_climo = add_coords(compute_climatology(load_var('vo')))
    for mon in [3, 6, 9, 12]:
        print(f'  month {mon:02d}')
        u_mon = add_valid_time(u_climo.sel(month=mon))
        v_mon = add_valid_time(v_climo.sel(month=mon))
        for seg in segments:
            print('    ' + seg.segstr)
            regridded = seg.regrid_velocity(
                u_mon, v_mon,
                regrid_suffix='spear_velocity',
                write=False
            )
            seg.to_netcdf(regridded, f'spear_uv_i{mon:02d}_climo')

