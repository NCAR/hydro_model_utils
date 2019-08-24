#  #!/glade/u/home/gangrade/miniconda3/envs/analysis

import xarray as xr
import glob, os

ncdir = '../output/'
finalfile = 'cali_huc12.merged.nc'

outfilelist = glob.glob((ncdir+'/*{}*.nc').format('output.CALI_huc12*')) # file list
outfilelist.sort()

nclist = []  #List for storing the converted SUMMA output files

for count, value in enumerate(outfilelist):
    ncconvert = xr.open_dataset(outfilelist[count])                   # Import netCDF file
    runoffdata = ncconvert['scalarTotalRunoff_mean'].values           # Extract scalarTotalRunoff_mean values
    runoffarray = xr.DataArray(runoffdata, dims=['time','hru'])       # Create array of scalarTotalRunoff_mean with 2 dimensions
    ncconvert = ncconvert.drop('scalarTotalRunoff_mean')              # Drop the original scalarTotalRunoff_mean variable
    ncconvert['scalarTotalRunoff_mean'] = runoffarray                 # Add the new array to original netCDF
    ncconvert['scalarTotalRunoff_mean'].attrs['long_name'] = "total runoff (mean)"
    ncconvert['scalarTotalRunoff_mean'].attrs['units'] = 'm s-1'
    nclist.append(ncconvert)

#Part 2: Converted all output files
ncconcat = xr.concat(nclist, dim='hru')

#Remove NANs because mizuRoute doesn't like it
FILL_VALUE = -9999.0
#ncconcat['pptrate'].encoding = {'_FillValue':FILL_VALUE }
ncconcat['time'].encoding = {'_FillValue':FILL_VALUE }
ncconcat['scalarTotalRunoff_mean'].encoding = {'_FillValue':FILL_VALUE }

ncconcat.to_netcdf(os.path.join(ncdir, finalfile), 'w')


