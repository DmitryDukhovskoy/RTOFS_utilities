import datetime as dt
import numpy as np
import pandas as pd
from pathlib import Path
import xarray

import sys
sys.path.append('../../setup/boundary')
from boundary import Segment

from write_spear_climatology_boundary import add_coords


def interpolate_segment_z(ds_in, varname, zname, depth_from, depth_to):
    ds = ds_in.copy()
    # Add depth_from to input data and set as main z coord
    ds['depth'] = ((zname, ), depth_from)
    ds = ds.swap_dims({zname: 'depth'})
    # Interpolate input data onto depth_to depth values
    interped = ds[varname].interp(depth=depth_to).ffill('depth') # need to ffill because glorys is shallow
    # Go back to zname as the main z coordinate
    interped[zname] = (('depth', ), np.arange(len(interped.depth)))
    interped = interped.swap_dims({'depth': zname})
    return interped


def load_glorys_climo(var, seg, glorys_var=None):
    if glorys_var is None:
        glorys_var = var
    glorys_climo = xarray.open_dataset(f'/work/acr/mom6/nwa12/glorys/updated/climatology/glorys_{glorys_var}_climo_{seg.num:03d}.nc')
    # Load GLORYS climo
    if var != glorys_var:
        glorys_climo = glorys_climo.rename({f'{glorys_var}_segment_{seg.num:03d}': f'{var}_segment_{seg.num:03d}'})
    # Add time information to match forecast climatology
    glorys_climo['time'] = (('month', ), pd.date_range('1993-01-01', periods=12, freq='1MS'))
    glorys_climo = glorys_climo.swap_dims({'month': 'time'})
    return glorys_climo


def load_spear_climo(var, seg, mstart): # note global mstart
    if var in ['uo', 'vo']:
        spear_climo = xarray.open_dataset(f'/work/acr/spear/climatology/spear_uv_i{mstart:02d}_climo_{seg.num:03d}.nc')
    else:
        spear_climo = xarray.open_dataset(f'/work/acr/spear/climatology/spear_{var}_i{mstart:02d}_climo_{seg.num:03d}.nc')
    spear_climo = spear_climo.drop('lead', errors='ignore')
    return spear_climo


def rep_year(ds):
    rep = ds.copy()
    rep['time'] = ds.time.to_pandas() + pd.DateOffset(years=1)
    ds = xarray.concat((ds, rep), dim='time')
    return ds


def main(ystart, mstart, ens):
    indir = Path('/work/acr/spear/processed/')
    outdir = Path(f'/work/acr/mom6/nwa12/forecast_input_data/boundary/{ystart}-{mstart:02d}')
    outdir.mkdir(exist_ok=True)
    hgrid = xarray.open_dataset('../../setup/grid/ocean_hgrid.nc')
    segments = [
        Segment(1, 'south', hgrid, output_dir=outdir),
        Segment(2, 'north', hgrid, output_dir=outdir),
        Segment(3, 'east', hgrid, output_dir=outdir)
    ]
    glorys_z = xarray.open_dataset('/work/acr/glorys/GLOBAL_MULTIYEAR_PHY_001_030/depths.nc').depth
    spear_z = xarray.open_dataset('/work/acr/spear/depths.nc').z_l

    print('Velocity')
    glorys_var = 'uv'
    u_forecasts = xarray.open_dataset(indir / f'ocean_z.monthly_mean.ens_{ens:02d}.uo.nc')
    v_forecasts = xarray.open_dataset(indir / f'ocean_z.monthly_mean.ens_{ens:02d}.vo.nc')
    for seg in segments:
        print(f'  {seg.num}')
        varstr = '{var}_segment_{seg.num:03d}'
        zstr = f'nz_segment_{seg.num:03d}'
        dzstr = 'dz_{var}_segment_{seg.num:03d}'

        glorys_climo = load_glorys_climo('uv', seg)
        # Interpolate glorys climatology onto SPEAR vertical grid
        interped = xarray.merge((
            interpolate_segment_z(glorys_climo, varstr.format(var='u', seg=seg), zstr, glorys_z.data, spear_z.data),
            interpolate_segment_z(glorys_climo, varstr.format(var='v', seg=seg), zstr, glorys_z.data, spear_z.data)
        ))
        # Duplicate for a second year, so that subtracting a forecast climatology that 
        # spans multiple years works
        interped = rep_year(interped)
        # Load SPEAR climatology for forecasts started in {mstart}
        spear_climo = load_spear_climo('uv', seg, mstart)
        # Pick the forecast for {ystart}, {mstart} and {ens}
        this_forecasts = []
        for forecasts in [u_forecasts, v_forecasts]:
            this_forecast = forecasts.sel(start=f'{ystart}-{mstart}').squeeze()
            this_forecast = add_coords(this_forecast)
            if 'z_l' in this_forecast.coords:
                this_forecast = this_forecast.rename({'z_l': 'z'})
            this_forecast = this_forecast.rename({'valid_time': 'time'}).swap_dims({'lead': 'time'}).drop('lead', errors='ignore')
            this_forecasts.append(this_forecast)
        
        usource = this_forecasts[0]
        vsource = this_forecasts[1]

        # Regrid to model boundary
        regridded = seg.regrid_velocity(
            usource, vsource,
            uvar='uo', 
            vvar='vo',
            regrid_suffix='spear_velocity',
            write=False
        )

        corrected = regridded.copy()
        # corrected = forecast - bias = forecast - (forecast_climo - glorys_climo)
        for v in ['u', 'v']:
            var = varstr.format(var=v, seg=seg)
            bias = spear_climo[var] - interped[var].sel(time=spear_climo.time)
            bias['time'] = regridded['time'] # assumption that months match (years differ)
            corrected[var] = corrected[var] - bias

        # Load the last glorys daily value before the forecast start
        # [this does not have logic to handle a start on January 1]
        last_glorys = (
            xarray.open_dataset(f'/work/acr/mom6/nwa12/glorys/updated/{glorys_var}_{seg.num:03d}_{ystart}.nc')
            .sel(time=dt.datetime(ystart, mstart, 1) - dt.timedelta(hours=12)) # GLORYS time is at noon
        )
        head = xarray.merge((
            interpolate_segment_z(last_glorys, varstr.format(var='u', seg=seg), zstr, glorys_z.data, spear_z.data),
            interpolate_segment_z(last_glorys, varstr.format(var='v', seg=seg), zstr, glorys_z.data, spear_z.data)
        ))
        # Find the time in the middle of the month after the end of the forecast
        # [this assumes that forecasts are 1 year long, e.g. from Jan1--Dec31]
        next_mid = regridded['time'][0].values + pd.DateOffset(years=1)
        # Find GLORYS climatology for the first forecast month
        # and duplicate it as the climatology of the month
        # after the end of the forecast (same month for 1 year forecasts)
        tail = interped.sel(time=spear_climo.time).isel(time=0)
        tail['time'] = next_mid

        for v in ['u', 'v']:
            dz = dzstr.format(var=v, seg=seg)
            head[dz] = corrected[dz].isel(time=0)
            tail[dz] = corrected[dz].isel(time=0)

        # Combine the forecast with the two padding values
        combined = xarray.concat((head, corrected, tail), dim='time')
        combined = combined.transpose(*(['time'] + [d for d in combined[varstr.format(var='u', seg=seg)].dims if d != 'time']))
        # Make sure that the calendar becomes gregorian in the output
        if 'calendar' in combined['time'].attrs:
            del combined['time'].attrs['calendar']        
        seg.to_netcdf(combined, f'forecast_v_i{ystart}-{mstart:02d}-e{ens:02d}')

    for var in ['ssh', 'so', 'thetao']:
        print(var)
        fname = f'ocean.monthly_mean.ens_{ens:02d}.{var}.nc' if var == 'ssh' else f'ocean_z.monthly_mean.ens_{ens:02d}.{var}.nc'
        all_forecasts = xarray.open_dataset(indir / fname)
        if var == 'ssh':
            glorys_var = 'zos'
        else:
            glorys_var = var
        for seg in segments:
            print(f'  {seg.num}')
            varstr = f'{var}_segment_{seg.num:03d}'
            zstr = f'nz_segment_{seg.num:03d}'
            dzstr = f'dz_{var}_segment_{seg.num:03d}'
            glorys_climo = load_glorys_climo(var, seg, glorys_var=glorys_var)
            # Interpolate glorys climatology onto SPEAR vertical grid
            if var == 'ssh':
                interped = glorys_climo[varstr].copy()
            else:
                interped = interpolate_segment_z(glorys_climo, varstr, zstr, glorys_z.data, spear_z.data)
            # Duplicate for a second year, so that subtracting a forecast climatology that 
            # spans multiple years works
            interped = rep_year(interped)
            # Load SPEAR climatology for forecasts started in {mstart}
            spear_climo = load_spear_climo(var, seg, mstart)
            # Pick the forecast for {ystart}, {mstart} and {ens}
            this_forecast = all_forecasts.sel(start=f'{ystart}-{mstart}').squeeze()
            # Add lat/lon to help regrid
            this_forecast = add_coords(this_forecast)
            if 'z_l' in this_forecast.coords:
                this_forecast = this_forecast.rename({'z_l': 'z'})
            this_forecast = this_forecast.rename({'valid_time': 'time'}).swap_dims({'lead': 'time'}).drop('lead', errors='ignore')
            # Regrid to model boundary
            regridded = seg.regrid_tracer(
                this_forecast,
                source_var=var,
                regrid_suffix='spear_tracer',
                write=False
            )
            # corrected = forecast - bias = forecast - (forecast_climo - glorys_climo)
            bias = spear_climo[varstr] - interped.sel(time=spear_climo.time)
            bias['time'] = regridded['time'] # assumption that months match (years differ)
            corrected = regridded.copy()
            corrected[varstr] = corrected[varstr] - bias
            # Load the last glorys daily value before the forecast start
            # [this does not have logic to handle a start on January 1]
            last_glorys = (
                xarray.open_dataset(f'/work/acr/mom6/nwa12/glorys/updated/{glorys_var}_{seg.num:03d}_{ystart}.nc')
                .sel(time=dt.datetime(ystart, mstart, 1) - dt.timedelta(hours=12)) # GLORYS time is at noon
            )
            if var != glorys_var:
                last_glorys = last_glorys.rename({
                    f'{glorys_var}_segment_{seg.num:03d}': f'{var}_segment_{seg.num:03d}',
                    # dzstr.replace(var, glorys_var): dzstr # this would be needed for 3D var with different name
                })
            if var == 'ssh':
                head = last_glorys[varstr]
            else:
                head = interpolate_segment_z(last_glorys, varstr, zstr, glorys_z.data, spear_z.data)
            head = head.to_dataset()
            if dzstr in corrected: # SSH does not have this
                head[dzstr] = corrected[dzstr].isel(time=0)
            # Find the time in the middle of the month after the end of the forecast
            # [this assumes that forecasts are 1 year long, e.g. from Jan1--Dec31]
            next_mid = regridded['time'][0].values + pd.DateOffset(years=1)
            # Find GLORYS climatology for the first forecast month
            # and duplicate it as the climatology of the month
            # after the end of the forecast (same month for 1 year forecasts)
            tail = interped.sel(time=spear_climo.time).isel(time=0)
            tail['time'] = next_mid
            tail = tail.to_dataset()
            if dzstr in corrected: # SSH does not have this
                tail[dzstr] = corrected[dzstr].isel(time=0)
            # Combine the forecast with the two padding values
            combined = xarray.concat((head, corrected, tail), dim='time')
            combined = combined.transpose(*(['time'] + [d for d in combined[varstr].dims if d != 'time']))
            # Make sure that the calendar becomes gregorian in the output
            if 'calendar' in combined['time'].attrs:
                del combined['time'].attrs['calendar']        
            seg.to_netcdf(combined, f'forecast_{var}_i{ystart}-{mstart:02d}-e{ens:02d}')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('YEAR', type=int)
    parser.add_argument('MONTH', type=int)
    parser.add_argument('ENS', type=int)
    args = parser.parse_args()    
    main(args.YEAR, args.MONTH, args.ENS)

