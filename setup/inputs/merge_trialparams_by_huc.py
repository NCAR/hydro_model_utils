# requires:  conda install -c conda-forge xarray dask distributed netcdf4

from dask.distributed import Client
import xarray as xr

#client = Client()  # it registers some stuff behind the scenes. you can print it (`print(client)`), 
                   # it will tell you what computational resources it has available. The point of 
                   # this line is to enable parallel reads.

                   # don't use if on a cheyenne interactive node!

in_template = '/glade/u/home/gangrade/WORK/Summa_Preprocessing/02_Reproject/ACT_HUC12.ParamTrial.label1.*.nc'
out_template = '/glade/u/home/gangrade/WORK/Summa_Preprocessing/02_Reproject/ACT_HUC12.ParamTrial.nc'
# years = range(2010, 2015)

ds = xr.open_mfdataset(in_template, concat_dim='hru')
ds.load().to_netcdf(out_template)
#ds.load().to_netcdf(out_template, unlimited_dims=['hru'])



#in_template = '/glade2/scratch2/andywood/SHARP/wreg/forcing/huc12_3hr/huc12_forcing_wUS.{0:4d}.3hr.*.nc'
#out_template = '/glade/u/home/jhamman/workdir/summa_junk/huc12_forcing_wUS.{0:4d}.3hr.nc'
#for y in years:
#    ds = xr.open_mfdataset(in_template.format(y), concat_dim='hru')
#    ds.load().to_netcdf(out_template.format(y), unlimited_dims=['time'])


