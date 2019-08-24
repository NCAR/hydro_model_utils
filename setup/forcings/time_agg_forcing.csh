#!/bin/csh
# aggregate forcings in time (1 hr to 3 hr)
# A. Wood

set InDir = ./nldas/
set OutDir = ./nldas_3hr/

module load cdo

foreach F (`\ls $InDir/forc*1980*`)

  echo "-----------"
  echo processing $F
  cdo -timselmean,3,0,0 -shifttime,-3600seconds $F temp.nc
  ncks -C -x -v time_bnds temp.nc $OutDir/$F:t:r.3hr.nc   # drop new variable, not needed
  ncks -A data_step.3hr.nc $OutDir/$F:t:r.3hr.nc  # update data_step
  \rm temp.nc

end

