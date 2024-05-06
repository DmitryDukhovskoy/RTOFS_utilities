"""
  Create nudging inverse time scale
  usage: write_damping_file.py -c config_nep.yaml -r N
  N = int days
"""
import xarray

def create_damping(ocean_static, rate):
    t_ds = ocean_static['wet'].copy()
    t_ds.name = 'Idamp'
    t_ds *= rate
    t_ds = t_ds.to_dataset()
    t_ds['Idamp'].attrs['units'] = 's-1'
    t_ds['Idamp'].attrs['cell_methods'] = 'time: point'
    return t_ds

if __name__ == '__main__':
    import argparse
    from pathlib import Path
    from yaml import safe_load
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config')
    parser.add_argument('-r', '--rate', type=int) # Note: int days only
    args = parser.parse_args()
    with open(args.config, 'r') as file: 
        config = safe_load(file)
    static = xarray.open_dataset(config['filesystem']['ocean_static'])
    model_input = Path(config['filesystem']['model_input_data']) / 'nudging'
    model_input.mkdir(exist_ok=True)
    damping = create_damping(static, 1 / (args.rate * 24 * 3600))
    encoding = {'Idamp': {'_FillValue': None}} 

    flout = f'damping_full_t_{args.rate:02d}.nc'
    print(f'Saving {model_input}/{flout}')
    damping.to_netcdf(
        model_input / flout,
        format='NETCDF3_64BIT',
        engine='netcdf4',
        encoding=encoding
    )

