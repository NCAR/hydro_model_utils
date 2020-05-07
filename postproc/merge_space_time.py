# script to combine split domain results (with 2 dimensions, hru and/or time)
# input details hardwired for now
# AWW 2019

# requires:  conda install -c conda-forge xarray dask distributed netcdf4
# or use built-in env, eg on cheyenne
# >module load python/3.6.8
# >source /glade/u/apps/opt/ncar_pylib/ncar_pylib.csh

# example input/output names
# inFnameTemplate = 'testout/3hr/camels_v2_3hr_3lyr_output_*_day.nc'
# outFname        = 'testout/3hr/camels_v2_3hr_3lyr_output_day.nc'
# -----------------------------------------------

import sys
#from dask.distributed import Client  (only if using paralellism)

# --------- arguments -----------
print("found %d args" % len(sys.argv))
if len(sys.argv) == 4:
    inFnameTemplate = sys.argv[1]
    outFname        = sys.argv[2]
    mergeDim        = sys.argv[3]  # e.g., 'hru'
    print("input %s" % inFnameTemplate)
    print("output %s" % outFname)

else:
    print("USAGE: %s <input_filename_template> <output_filename> <mergeDim>" % sys.argv[0])
    sys.exit(0)

# now import slower packages 
import xarray as xr

# --------- code -----------

# identify dimension names in files
unlimDim = 'time'

# the next line enables parallel reads, but it's fast without them
# client = Client() # this registers some stuff behind the scenes. you 
                    # can print it (`print(client)`) to see what comp. resources 
                    # are available. 

ds = xr.open_mfdataset(inFnameTemplate, concat_dim=mergeDim)
ds.load().to_netcdf(outFname, unlimited_dims=[unlimDim])

# DONE





## ===== not used here but the following allows looping over years if needed

#years = range(2007, 2020)  # NOTE: last year not included

#in_template  = 'data/gmet_huc12/tmpsplit/ens_forc.WEST.huc12.{0:4d}.001.*.nc'
#out_template = 'data/gmet_huc12/ens_forc.WEST.huc12.{0:4d}.001.nc'

#for y in years:
#    print("working on year %d" % y)
#    ds = xr.open_mfdataset(in_template.format(y), concat_dim='hru')
#    ds.load().to_netcdf(out_template.format(y), unlimited_dims=['time'])

