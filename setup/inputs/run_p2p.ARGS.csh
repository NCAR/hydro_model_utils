#!/bin/csh
# A. Wood, 2018
# run python script to apply poly2poly to nldas vector data:  attributes

set start_poly = $1 # first in file is 0
set end_poly   = $2 # maximum of number of target polygons - 1
set label      = $3 # sets apart parallel output files for different poly ranges
set modelName  = R18_HUC12
# set Mapping    = /glade/u/home/andywood/proj/SHARP/wreg/p2p/mapping.nldas_hru_id.to.ReclDomainTrim.v0.nc
# set IDMapFile  = /glade/u/home/andywood/proj/SHARP/wreg/p2p/IDMap.mc_nldas_forcing.hruId_to_index.txt

# set Mapping    = /glade/u/home/gangrade/WORK/Summa_Preprocessing/01_Mapping/mapping_nldas_to_HUC12_ACT.nc
# set Mapping    = /glade/u/home/gangrade/WORK/Summa_Preprocessing/01_Mapping/mapping_nldas_to_milk_stmary.nc
# set Mapping    = /glade/u/home/gangrade/WORK/Summa_Preprocessing/01_Mapping/mapping_nldas_to_WESTHUC12.nc
# set Mapping    = /glade/u/home/gangrade/WORK/Summa_Preprocessing/01_Mapping/mapping_nldas_to_3150101.nc
set Mapping    = /glade/u/home/gangrade/WORK/Summa_Preprocessing/01_Mapping/mapping_nldas_to_R18_HUC12.nc

set IDMapFile  = /glade/u/home/andywood/proj/SHARP/wreg/p2p/IDMap.new_nldasID_to_index.filled.txt


set HruType    = int     # actually 'i4' or 'int64'
set Python     = python2.7

# set InFl       = /glade/u/home/andywood/proj/SHARP/wreg/nldas_config/summa_run/settings/nldasAttributes_vector.nc
set InFl       = /glade/u/home/gangrade/WORK/LIS_summa/nldasAttributes_vector.v1.test.nc
set OutFl      = /glade/u/home/gangrade/WORK/Summa_Preprocessing/02_Reproject/$modelName.Attributes_HRU.$label.nc
#set OutFl      = ./wUShuc12_Attributes.$label.nc
set OutFl2     = /glade/u/home/gangrade/WORK/Summa_Preprocessing/02_Reproject/$modelName.Attributes_GRU.$label.nc

# choose one of the following:

# generate attributes with 'hru' index (~40 minutes)
$Python poly2poly.map_vector_attributes.py $Mapping $InFl $OutFl $HruType $IDMapFile $start_poly $end_poly 

# generate additional variable with 'gru' index (to be added later)  (~5 minutes)
# $Python poly2poly.map_vector_GRUvars.py $Mapping $InFl $OutFl2 $HruType $IDMapFile $start_poly $end_poly
