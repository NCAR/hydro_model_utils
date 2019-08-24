#!/bin/csh
# A. Wood
# make monthly forcing files into yearly ones

set Y = 2010
while ($Y <= 2014)
  echo processing $Y
  ncrcat nldas_3hr/*$Y*.3hr.nc nldas_3hr_yrly/nldasForcing_$Y.3hr.nc
  @ Y ++
end


#set Y = 1999
#  ncrcat nldas_3hr/*$Y*.3hr.nc nldas_3hr_yrly/nldasForcing_$Y.3hr.nc
#exit


