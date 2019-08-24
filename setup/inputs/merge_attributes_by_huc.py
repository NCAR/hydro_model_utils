# requires:  conda install -c conda-forge xarray dask distributed netcdf4
#  run with 'mypython2' (miniconda)

from dask.distributed import Client
import xarray as xr

#client = Client()  # it registers some stuff behind the scenes. you can print it (`print(client)`), 
                   # it will tell you what computational resources it has available. The point of 
                   # this line is to enable parallel reads.
                   # don't use if on a cheyenne interactive node!


in_template = '/glade/work/gangrade/Summa_Preprocessing/02_Reproject/WESTHUC12.Attributes_HRU.*.nc'
out_template = '/glade/work/gangrade/Summa_Preprocessing/02_Reproject/WESTHUC12.Attributes_HRU.all.nc'

ds = xr.open_mfdataset(in_template, concat_dim='hru')
ds.load().to_netcdf(out_template)
