#!/bin/csh -x
# merge split domain summa run into usable output file
# need to update filenames for each run
# A. Wood, mar 2020
# to run on cheyenne/casper, first
#   module load python/3.6.8
#   source /glade/u/apps/opt/ncar_pylib/ncar_pylib.csh
#   module load nco

# --- settings 
set dataDir  = ../output/nldas_test/
set split    = $dataDir/nldas_v2_test_\*day.nc
set split_g  = $dataDir/nldas_v2_test_\*day.nc.grus   # if needed
set merged   = $dataDir/nldas_v2_test.day.nc
set merged_g = $dataDir/nldas_v2_test.day.GRU.nc     # if needed
set runoff   = $dataDir/nldas_v2_test.day.runoff.nc
set gruVars  = averageRoutedRunoff_mean  # "none" if no gru vars are in split output; else a comma-sep list

# --- if grus are present in split files, need to separate out those variables, otherwise merge existing files
if ( gruVars == "none") then

  # --- merge all split files into one by running python script
  echo merging hrus
  python /glade/u/home/andywood/proj/aist/nldas/summa/camels/postproc/merge_space_time.py "$split" "$merged" "hru"

  # --- extract runoff file for routing (if needed)
  echo extracting runoff
  ncks -h -v scalarTotalRunoff_mean $merged $runoff

else

  # --- separate hru & gru variables (before merging)
  foreach F ($split)
    ncks -h -v $gruVars $F $F.grus
  end
  foreach F ($split)
    ncks -h -x -v $gruVars $F $F.hrus
    mv $F.hrus $F
  end

  # --- merge all split files into hru & gru combined files using python script
  echo merging hrus ... 
  python /glade/u/home/andywood/proj/aist/nldas/summa/camels/postproc/merge_space_time.py "$split" "$merged" "hru" 
  echo merging grus ... 
  python /glade/u/home/andywood/proj/aist/nldas/summa/camels/postproc/merge_space_time.py "$split_g" "$merged_g" "gru"

  # --- extract runoff file for routing (if needed)
  echo extracting runoff
  ncks -h -v scalarTotalRunoff_mean $merged $runoff

  # --- for gru-file, rename dim to 'hru' and add back into main hru file (works if hrus=grus)
  echo merging gru vars back into main hru file
  ncrename -d gru,hru $merged_g $dataDir/tmp.nc
  ncks -h -A $dataDir/tmp.nc $merged
  \rm $dataDir/tmp.nc

endif

# --- clean up run parts (uncomment f confident)
#\rm $split*

