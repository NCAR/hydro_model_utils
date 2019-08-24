#!/bin/csh
# A. Wood, 2018
# make a job list to run 36
# set for multiple years and months

set nHru = 2863   # specific to each domain
set nCoresPerNode = 36   # specific to cluster
set jobsPerCore = `echo $nHru $nCoresPerNode | awk '{printf("%d",$1/$2+1)}' `
set jobList = joblist.run_summa.txt
echo jobs per core = $jobsPerCore

echo -n > $jobList
set N = 0
while ($N < $nCoresPerNode)
  set N1 = `echo $N | awk '{print $1*'$jobsPerCore'+1}' `
  set N2 = `echo $N1 $jobsPerCore $nHru | awk '{if(($1+$2)>$3){print $3-$1+1}else{print $2}}' `
  echo ./summa.exe -g $N1 $N2 -r never -m fileManager.txt >> $jobList
  @ N ++
end
