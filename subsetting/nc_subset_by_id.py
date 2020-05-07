#!/usr/bin/env python
# Script to subset a netcdf file based on a list of 
# IDs

import argparse
from datetime import datetime
import os
import sys
import xarray as xr

def process_command_line():
    '''Parse the commandline'''
    parser = argparse.ArgumentParser(description='Script to subset a netcdf file based on a list of IDs.')
    parser.add_argument('ncfile',
                        help='path of netcdf file that will be subset.')
    parser.add_argument('idfile', 
                        help='path of file with list of ids.')
    parser.add_argument('opath',
                        help='directory where subsetted file will be written.')
    args = parser.parse_args()
    return(args)


# main
if __name__ == '__main__':
    # process command line
    args = process_command_line()

    # read the IDs to subset
    with open(args.idfile) as f:
        ids = [int(x) for x in f]

    # subset dimension
    subset_dim = 'hru'

    # ingest the netcdf file
    ds = xr.open_dataset(args.ncfile)

    # determine which variables do not vary along any of the subsetted dimensions
    vars_without = []
    for var in ds.variables:
        if not subset_dim in ds[var].dims:
            vars_without.append(var)

    # take the other variables
    vars_with = list(set(ds.variables).difference(set(vars_without)))

    # subset the netcdf file based on the hruId
    ds_with = ds[vars_with]
    ds_subset = ds_with.where(ds_with.hruId.isin(ids), drop=True)

    # merge the with and without
    ds_subset = ds_subset.merge(ds[vars_without])

    # make sure that the subsetted types are the same as the original ones
    for var in ds_subset.variables:
        ds_subset[var] = ds_subset[var].astype(ds[var].dtype)

    # update the history attribute
    history = '{}: {}\n'.format(datetime.now().strftime('%c'), 
                                ' '.join(sys.argv))
    if 'history' in ds_subset.attrs:
        ds_subset.attrs['history'] = history + ds_subset.attrs['history']
    else:
        ds_subset.attrs['history'] = history

    # Write to file
    ofile = os.path.join(args.opath, os.path.basename(args.ncfile))
    ds_subset.to_netcdf(ofile)

    # Write IDs from the ID file that were not in the NetCDF file to stdout
    missing = set(ids).difference(set(ds_subset.hruId.values))
    if missing:
        print("Missing IDs: ")
        for x in missing:
            print(x)
