#!/bin/csh -e
# S. Gangrade, 2019 
# makes job lists and automatically submit jobs to create forcings starting from year 1980 to 2012

# \rm joblist.run_p2p.txt

set run = 1;
set Yr = 1980
while ($Yr <= 2015)
  
  set N = 0
  while ($N < 36)
    set mm = 1
    while ($mm <= 12)
      # month = printf("%03d",$mm)

      set month = `echo $Yr $mm | awk '{printf("%04d",$1)}''{printf("%02d",$2)}'`
      # set xx = "$(printf "%02d" 5)"
      # echo $foo
      echo $month
      echo /glade/work/gangrade/Summa_Preprocessing/02_Reproject/forcings/run_p2p.ARGS.csh $month 0 4452 >> joblist.run_p2p.$run.txt
      @ mm ++
      @ N ++
    end
	@ Yr ++
  end
  @ run ++
  # cd /glade/u/home/gangrade/WORK/Summa_Preprocessing/02_Reproject/forcings 
  # qsub lsf.run_p2p.SG.pbs
  # \rm joblist.run_p2p.txt
end
