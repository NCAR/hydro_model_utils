#!/bin/csh
# A. Wood, 2018
# run python script to apply poly2poly to forcing data


set year       = 2010
set start_poly = 0
set end_poly   = 1
set label      = 999

set Mapping    = /glade/u/home/andywood/proj/SHARP/wreg/p2p/mapping.nldas_hru_id.to.ReclDomainTrim.v0.nc
set IDMapFile  = /glade/u/home/andywood/proj/SHARP/wreg/p2p/IDMap.mc_nldas_forcing.hruId_to_index.txt
set HruType    = int     # actually 'i4'
set InFl       = /glade/u/home/andywood/scratch/SHARP/wreg/forcing/nldas_3hr_yrly/nldasForcing_$year.3hr.nc
#set OutFl      = /glade/u/home/andywood/scratch/SHARP/wreg/forcing/huc12_3hr/huc12_forcing_wUS.$year.3hr.$label.nc
set OutFl      = ./testout.nc

set Python     = /glade/u/home/andywood/mylibs/miniconda2/bin/python
set PyScript   = /glade/u/home/andywood/proj/SHARP/wreg/p2p/forcings/poly2poly.map_vector_timeseries.py


$Python $PyScript $Mapping $InFl $OutFl $HruType $IDMapFile $start_poly $end_poly
