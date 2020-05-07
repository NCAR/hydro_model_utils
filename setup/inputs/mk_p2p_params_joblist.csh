#!/bin/csh
# run python script to apply poly2poly to parameter data

# before running jobs need to load python env
#  module load python/2.7.16; ncar_pylib

set nOutPolys = 54741
set nCores    = 36
set jobsPerCore = `echo $nOutPolys $nCores | awk '{printf("%d",$1/$2+1)}' `
set jobList     = ./joblist.p2p_trialparams.txt
echo jobs per core = $jobsPerCore


set Mapping   = maps/mapping.NLDAS125_to_wUShuc12x.nc
set InFl      = ../test_params//nldasParamTrial.nc
set OutFl     = ./huc_params/wUShuc12_ParamTrial.$label.nc
set HruType   = int  # actually 'i4'
set IDMapFile = IDtoIndex.nldas_conus.mapping.txt
set Python    = /Users/andywood/misc/miniconda2/bin/python2.7

echo -n > $jobList
set N = 0
while ($N < $nCores)
  set N1 = `echo $N | awk '{print $1*'$jobsPerCore'+1}' `
  set N2 = `echo $N1 $jobsPerCore $nHru | awk '{if(($1+$2)>$3){print $3-$1+1}else{print $2}}' `
  echo $Python poly2poly.map_vector_trialparams.py $Mapping $InFl $OutFl $HruType $IDMapFile $N1 $N2 $start_poly $end_poly $varname

  @ N ++
end



#$Python poly2poly.map_vector_trialparams.py $Mapping $InFl $OutFl $HruType $IDMapFile $start_poly $end_poly $varname

