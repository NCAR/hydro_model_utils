#!/bin/csh
# run python script to apply poly2poly to attribute data

set start_poly = 0  # starts at 0
set end_poly = 1
set label = 001

set Mapping   = ../p2p_maps/mapping.nldas_hru_id.to.ReclDomainTrim.v0.nc
set InFl      = ../test_params/nldasAttributes_vector.nc
set OutFl     = ./huc_params/wUShuc12_Attributes.$label.nc
set HruType   = int  # actually 'i4'
set IDMapFile = IDtoIndex.nldas_conus.mapping.txt
set Python    = /Users/andywood/misc/miniconda2/bin/python2.7


echo $Python poly2poly.map_vector_attributes.py $Mapping $InFl $OutFl $HruType $IDMapFile $start_poly $end_poly

$Python poly2poly.map_vector_attributes.py $Mapping $InFl $OutFl $HruType $IDMapFile $start_poly $end_poly


set OutFl     = ./huc_params/wUShuc12_GRUvar.$label.nc
$Python poly2poly.map_vector_GRUvars.py $Mapping $InFl $OutFl $HruType $IDMapFile $start_poly $end_poly
