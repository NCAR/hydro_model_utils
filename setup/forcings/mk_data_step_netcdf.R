#!/usr/bin/env Rscript
# run these commands in R to create a single scalar variable to add to the forcing files
# A. Wood, 2018

library(ncdf4)

#data_step <- 10800 # 3 hrs
data_step <- 3600  # 1 hr
var.def   <- ncvar_def("data_step","seconds",dim=list(),missval=NULL,
                       "data step length in seconds",prec="double")
ncout     <- nc_create("data_step.1hr.nc",var.def)
ncvar_put(ncout,var.def,data_step)

nc_close(ncout)
