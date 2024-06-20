from pathlib import Path
import xarray


indir = Path('/work/acr/mom6/nwa12/glorys/updated')
outdir = Path('/work/acr/mom6/nwa12/glorys/updated/climatology')
outdir.mkdir(exist_ok=True)

for var in ['zos', 'uv', 'thetao', 'so']:
    print(var)
    for seg in [1, 2, 3]:
        print(f'  {seg}')
        source = xarray.open_dataset(indir / f'{var}_{seg:03d}.nc').sel(time=slice('1993', '2019'))
        monthly = source.groupby('time.month').mean('time')
        monthly.to_netcdf(outdir / f'glorys_{var}_climo_{seg:03d}.nc')


