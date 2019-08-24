#!/bin/csh
# A. Wood, 2018
# run python script to apply poly2poly to forcing data


set year       = $1
set start_poly = $2
set end_poly   = $3
set label      = R18_HUC12

# set Mapping    = /glade/u/home/andywood/proj/SHARP/wreg/p2p/mapping.nldas_hru_id.to.ReclDomainTrim.v0.nc
# set IDMapFile  = /glade/u/home/andywood/proj/SHARP/wreg/p2p/IDMap.mc_nldas_forcing.hruId_to_index.txt
# set Mapping    = /glade/u/home/gangrade/WORK/Summa_Preprocessing/01_Mapping/mapping_nldas_to_HUC12_ACT.nc
# set Mapping    = /glade/u/home/gangrade/WORK/Summa_Preprocessing/01_Mapping/mapping_nldas_to_milk_stmary.nc
# set Mapping    = /glade/u/home/gangrade/WORK/Summa_Preprocessing/01_Mapping/mapping_nldas_to_WESTHUC12.nc
# set Mapping    = /glade/u/home/gangrade/WORK/Summa_Preprocessing/01_Mapping/mapping_nldas_to_3150101.nc
set Mapping    = /glade/u/home/gangrade/WORK/Summa_Preprocessing/01_Mapping/mapping_nldas_to_R18_HUC12.nc


set IDMapFile  = /glade/u/home/andywood/proj/SHARP/wreg/p2p/IDMap.new_nldasID_to_index.filled.txt
set HruType    = int     # actually 'i4'
# set InFl       = /glade/u/home/ridwan/work/SUMMA/LIS_summa/lis_forcing/nldasForcing_$year.nc
set InFl       = /glade/u/home/andywood/rhap/common_data/nldas2_hr_summa/nldasForcing_$year.nc
set OutFl      = /glade/u/home/gangrade/WORK/Summa_Preprocessing/02_Reproject/$label.forcing.$year.nc

set Python     = python2.7
set PyScript   = /glade/u/home/gangrade/WORK/Summa_Preprocessing/02_Reproject/forcings/poly2poly.map_vector_timeseries.py


$Python $PyScript $Mapping $InFl $OutFl $HruType $IDMapFile $start_poly $end_poly
