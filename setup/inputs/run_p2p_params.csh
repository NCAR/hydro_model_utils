#!/bin/csh
# run python script to apply poly2poly to parameter data

set start_poly = 48850
set end_poly = 48852
set label = 001
set varname = 'frozenPrecipMultip'

set Mapping   = ../p2p_maps/mapping.nldas_hru_id.to.ReclDomainTrim.v0.nc
set InFl      = ../test_params//nldasParamTrial.nc
set OutFl     = ./huc_params/wUShuc12_ParamTrial.$label.nc
set HruType   = int  # actually 'i4'
set IDMapFile = IDtoIndex.nldas_conus.mapping.txt
set Python    = /Users/andywood/misc/miniconda2/bin/python2.7


echo $Python poly2poly.map_vector_trialparams.py $Mapping $InFl $OutFl $HruType $IDMapFile $start_poly $end_poly $varname
#$Python poly2poly.map_vector_trialparams.py $Mapping $InFl $OutFl $HruType $IDMapFile $start_poly $end_poly $varname

