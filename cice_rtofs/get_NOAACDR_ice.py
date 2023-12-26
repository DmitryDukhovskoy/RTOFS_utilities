"""
  Sea Ice Concentration NSIDC Near-Real-Time Climate Data Record V2, Arctic
  https://nsidc.org/data/user-resources/help-center/how-access-noaansidc-data-sets-polarwatch

 works from coaps machines but not from hera - firewall 
"""
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

def point_to_dataset(dataset_id, \
            base_url='https://polarwatch.noaa.gov/erddap/griddap'):
    base_url = base_url.rstrip('/')
    full_url = '/'.join([base_url, dataset_id])
    return nc.Dataset(full_url)
 

# 'nsidcG02202v4nhmday' is the unique ID of our interested data 
# from PolarWatch ERDDAP data server
da = point_to_dataset('nsidcG10016v2nh1day')

# prints variable names
da.variables.keys()

# set mapping crs to Cartopy's North Polar Stereo graphic
crs_epsg = ccrs.NorthPolarStereo(central_longitude=-45)

# set figure size
fig = plt.figure(figsize=[10, 10])

# set the map projection and associated boundaries
ax = plt.axes(projection = crs_epsg)
ax.set_extent([-3850000.0, 3750000.0, -5350000, 5850000.0],crs_epsg)
ax.coastlines()
ax.add_feature(cfeature.LAND)

# set the data crs using 'transform' 
# set the data crs as described in the netcdf metadata
cs = ax.pcolormesh(da['xgrid'], da['ygrid'], da['cdr_seaice_conc_monthly'][0][:] , 
                   cmap=plt.cm.Blues,  transform= ccrs.NorthPolarStereo(true_scale_latitude=70, central_longitude=-45)) #transform default is basemap specs

fig.colorbar(cs, ax=ax, location='bottom', shrink =0.8)
ax.set_title('Ice Concentration using Cartopy projection NorthPolarStereo()')

plt.show()

