#!/usr/bin/env python
''' Process timeseries grid into mean areal timeseries for arbitrary polygons
    Depends on mapping file from poly2poly.py and associated scripts by
    K. Sampson, NCAR '''
# Notes:  AWW: on yellowstone, needs 'module load netcdf4python'
#              on cheyenne, module load netcdf/2.7.16 and then type 'ncar_pylibs'
#                to set environment
#
# First Coder:    Naoki Mizukami, Feb 2014
# Modifications:
#   AWW Feb 2014, Documented code and adapted to work on BCSD5-NLDAS to
#                 HUC4 conversion
#   AWW Jan 2015, Adapted to ens. forc. generation
#                 writes out pcp, tmax, tmin, tmean, trange
#                 and write all out into one file (so has one fewer arg)
#          updated print statements to python 3 style, also raise statement
#   EC Jul 2016, Mod. to work with new hru wgtnc format, which
#          uses i_index and j_index instead of latitude and longitude
#   AWW Jul 2016, flipped j_index, shifted to match 0 range start, and altered
#          various loop indices to be more intuitive; improved messaging
#   EC Nov 2016, double-checked j_index flip and elaborated on related
#                comment for clarity. limited line lengths to 78 chars.
# =========================================================================

import sys
import os
import time
import getopt
import numpy as np
import netCDF4 as nc4

############################################
#              Class                       #
############################################
class wgtnc:
    """object of basin netCDF including hru and areal weight/i_index/j_index
       of others polygons  """
    def __init__(self,nc):
        """Initialization """
        self.nc=nc

    def getWgtHru(self,hru):
        """For given hru id, get weights of the intersected polygons and associated
           i_index, j_index"""
        wgtAll = getNetCDFData(self.nc, 'weight')
        hruAll = getNetCDFData(self.nc, 'IDmask')
         # EC get i_index --- LON direction
        i_indexAll = getNetCDFData(self.nc, 'i_index')
        # EC get j_index --- LAT direction
        j_indexAll = getNetCDFData(self.nc, 'j_index')
        hruList = self.getHruID() # EC-get hru id list
        # EC-get index corresponding to hru
        Index = [i for i, x in enumerate(hruAll) if x == hru]
        self.wgt = list(wgtAll[Index]) #Get wgt list for hru
        self.i_index = list(i_indexAll[Index]) #Get i_index list for hru
        self.j_index = list(j_indexAll[Index]) #Get j_index list for hru

        return (self.wgt, self.i_index, self.j_index)

    def getHruID(self):
        """ get hru ID list of basin"""
        # EC get polyid for hru list
        self.hruid = list(getNetCDFData(self.nc, 'polyid'))
        return self.hruid

########################################################################
#                                Modules                               #
########################################################################

def getNetCDFData(fn, varname):
    """Read <varname> variables from NetCDF <fn> """
    f = nc4.Dataset(fn,'r')
    data = f.variables[varname][:]
    f.close()
    return data

def getNetCDFDim(fn, varname):
    """ Get dimension name used in <varname> from NetCDF <fn> """
    f = nc4.Dataset(fn,'r')
    dim = f.variables[varname].dimensions

    data = dict()
    for dimname in dim:
        data[dimname]=f.variables[dimname][:]
    f.close()
    return data

# custom write function to handle ensemble forcing grid vars pcp, t_mean
# and t_range
def writeNetCDFData(fn, prcp, matavg, matrng, times, hrus, hru_type):
    """ Write <vars>[time,hru] array in netCDF4 file,<fn> and variable of
        <varname> """
    ncfile = nc4.Dataset(fn,'w',format='NETCDF4')

    # Data should be 2D [time x hru]
    dim2size = prcp.shape[1]
    # time axis - record dimension
    dim_1 = ncfile.createDimension('time', None)
    dim_2 = ncfile.createDimension('hru', dim2size)  # hru axis
    # EC add option for string or integer HRUs; AW for int64 and changed name to hruId
    if hru_type == 'str':
        # string HRU
        max_strlen = 12  # EC
        # EC: added string length for hrus
        dim_3 = ncfile.createDimension('strlen', max_strlen)
        hru = ncfile.createVariable('hruId', 'S1', ('hru', 'strlen')) 
    elif hru_type == 'int':
        # integer HRU
        hru = ncfile.createVariable('hruId', 'i4', ('hru', ))   
    elif hru_type == 'int64':
        # integer HRU
        hru = ncfile.createVariable('hruId', 'i8', ('hru', ))  
    else:
        print("hru type not recognized")
        sys.exit(1)

    # Define dimensions variable to hold the time var
    t = ncfile.createVariable('time', 'i4', ('time', ))

    # Define 2D variables to hold the data
    map = ncfile.createVariable('prcp','f4',('time','hru',))
    matx = ncfile.createVariable('tmax','f4',('time','hru',))
    matn = ncfile.createVariable('tmin','f4',('time','hru',))

    # Populate netcdf file variables
    t[:] = times   # NM added
    # EC/AW HRUs as integer or string or int64
    if hru_type == 'str':
        hru[:] = nc4.stringtochar(np.asarray(hrus,
                                  dtype='S{}'.format(max_strlen))) # EC
    elif hru_type == 'int':
        hru[:] = np.asarray(hrus, dtype='int')
    elif hru_type == 'int64':
        hru[:] = np.asarray(hrus, dtype='int64') # might need to try dtype=np.int64 (no 's)
    else:
        print("hru type not recognized")
        sys.exit(1)

    map[:,:]  = prcp   # for mean area precip (map), etc.
    matx[:,:] = (matavg + matrng/2)
    matn[:,:] = (matavg - matrng/2)

    # Attributes (should take more from input netcdf)
    # set time axis
    t.time_origin = '1970-JAN-01 00:00:00'
    t.title = 'Time'
    t.long_name = 'Time axis'
    t.units = 'seconds since 1970-01-01 00:00:00'
    t.calendar = 'Gregorian'
    t.axis = 'T'
    # set hru axis
    hru.title = 'elev zone ID'
    hru.long_name = 'elev zone axis'
    hru.axis = 'X'

    # set variable attributes
    map.associate = 'time hru'
    map.units = 'mm'
    map.setMissing = '1.e+20'
    map.axis = 'TX'

    matx.associate = 'time hru'
    matx.units = 'degrees C'
    matx.setMissing = '1.e+20'
    matx.axis = 'TX'

    matn.associate = 'time hru'
    matn.units = 'degrees C'
    matn.setMissing = '1.e+20'
    matn.axis = 'TX'

    # Write basic global attribute
    ncfile.history = 'Created ' + time.ctime(time.time())
    ncfile.source = os.path.dirname(os.path.abspath(__file__))+__file__[1:]

    ncfile.close()

def compAvgVal(nc_wgt, nc_in, varname):
    """Compute areal weighted avg value of <varname> in <nc_in> for each hru
       based on hru's weight in <nc_wgt>"""

    # instantiate current wgtnc object -- creates object
    wgt = wgtnc(nc_wgt)
    hruIDs = wgt.getHruID() # get hruID list

    dataVal = getNetCDFData(nc_in,varname)     # Get data value

    print("-------------------")
    print("INFO: var dataVal[0,0,0]: %s %f" % (varname, dataVal[0,0,0]))

    # AW: could speed this up by not importing this after the first time
    timeVal = getNetCDFData(nc_in,'time')

    dim1size = dataVal.shape[0]  # time dimension
    gridrows = dataVal.shape[1]  # grid rows dim (y)
    gridcols = dataVal.shape[2]  # grid cols dim (x)
    print("data shape",dataVal.shape)

    # AW: next two are not needed, but can be used during debugging (leave in)
    # GRID READ so [0,0] is LL corner
    #latVal = getNetCDFData(nc_in,'latitude')
    #lonVal  = getNetCDFData(nc_in,'longitude')

    # ==== LOOP over output polygons ====
    # !! this is inefficient - 
    # !!   weight/mapping info should be done once outside this call
    
    for h in range(start_poly,(end_poly+1)):        

    # Initialize data storage for weighted average time series for hru
    #for h in range(len(hruIDs)):
        # Get list of wgt, i_index, and j_index for corresponding hru
        (wgtval, i_index, j_index)=wgt.getWgtHru(hruIDs[h])

        ncells = len(wgtval)  # useful limit for several loops in this block

        # AW: need to (a) subtract 1 since meteorological grid indices in python
        #   start at [0,0] and input i/j indices start at [1, 1]; and
        #   (b) flip j-index because mapping file has northwest at [1, 1], so
        #   that's: j_index = gridrows - j_index

        for c in range(0,ncells):
            j_index[c] = gridrows - j_index[c]
            i_index[c] = i_index[c] - 1

        # ---- Calculate time series of weighted value for input polygon  ----

        # first, reset weights for the missing cells to zero
        sumweights = 0
        numvoid    = 0
        validcell  = 0
        # loop through cells(ie, c) assigned to this HRU
        for c in range(0,ncells):
            validcell = validcell + 1
            # EC: instead of lat -- y-index from UL corner (increments N->S)
            # instead of lon -- x-index from UL corner (increments W->E)
            yj = np.int(j_index[c])
            xi = np.int(i_index[c])
            # check if first timestep for cell is missing
            # if dataVal[1,yj,xi] > 1.e+19:
            # exclude cells with bad indexes (check why these occur in poly2poly
            if yj >= gridrows or xi >= gridcols or dataVal[1,yj,xi] > 1.e+19:
                wgtval[c] = 0
                numvoid = numvoid + 1
                print("WARNING: found bad value or bad indexes (i,j): %d %d" % (xi,yj))
            sumweights = sumweights + wgtval[c]

        print("-------------------")
        print(" var %s hru %d %s: num voids = %d out of %d" % (varname, h,
                                                               hruIDs[h],
                                                               numvoid,
                                                               validcell))
        print(" new sum of weights is %f" % sumweights)

        if sumweights > 0:
            # renormalize weights (won't work if all cells void)
            print(" normalizing...")
            for c in range(0,ncells):
                wgtval[c] = (wgtval[c]/sumweights)

            # now apply each wgt to matching cell data for all times
            flag = 0 # to track first time though
            for c in range(0,ncells):
                #print("working on cell c: %d" % (c))
                yj = np.int(j_index[c])
                xi = np.int(i_index[c])

                # weight values -- this is a vector with time
                Val1=wgtval[c]*dataVal[:,yj,xi].reshape(-1,1)

                if c == 0:
                    # initialize same sized vector
                    oriVal = Val1*0   # start with zero vector, won't affect sum

                # weight is above zero, add data from new cell
                if wgtval[c]>0:
                    oriVal=np.concatenate((oriVal, Val1), axis=1)

            Val2=oriVal.sum(axis=1).reshape(-1,1)   # sum up weighted HRU cell
            del oriVal

            # block for setting voids in HRUs with no good data
            if numvoid==validcell:
                print("ALL cells (%d of %d)" % (numvoid, validcell) +
                      " in hru %d are void, outputting voids" % (h))

                # assign voids (maybe not needed if using Val2 from above)
                for tm in range(len(Val2)):
                    Val2[tm] = 1.e+20

            # now store time vector of HRU results
            #if h == 0:
            #    wgtedVal = Val2
            #else:
            #    wgtedVal=np.concatenate((wgtedVal, Val2), axis=1)

            if 'wgtedVal' in locals():
                wgtedVal = np.concatenate((wgtedVal, Val2), axis=1) 
                activeOutPolyIDs = np.append(activeOutPolyIDs, hruIDs[h])
            else:
                # first time through
                wgtedVal = Val2
                activeOutPolyIDs = hruIDs[h]

        else:
            print("WARNING: no valid input contributing areas found")
            print("WARNING: all cells (%d of %d)" % (numvoid, validcell) +
                  " in hru %d are void or missing: omitting HRU from output" % (h))


    #return timeVal,hruIDs[start_poly:(end_poly+1)],wgtedVal
    return timeVal,activeOutPolyIDs,wgtedVal


############################################
#                Main                      #
############################################
use = '''
Usage: %s -[h] <weight_netCDF> <input_netCDF> <output_netCDF>
               <hru_type (int or str)> <start_poly:0- > <end_poly>
        -h  help
'''
if __name__ == '__main__':

    def usage():
        sys.stderr.write(use % sys.argv[0])
        sys.exit(1)
    try:
        (opts, args) = getopt.getopt(sys.argv[1:], 'h')
    except getopt.error:
        usage()

    verbose = False
    grid_info = False
    proj_info = True
    for (opt,val) in opts:
        if opt == '-h':
            usage()
        elif opt == '-v':
            verbose = True
        else:
            #raise OptionError, opt
            raise OptionError(opt)
            usage()

    if len(args) == 6:
        nc_wgt     = args[0]
        nc_in      = args[1]
        nc_out     = args[2]
        hru_type   = args[3]
        start_poly = int(args[4])  # starts at 0 (correspond w/ netcdf conventions)
        end_poly   = int(args[5])  # max possible value nOutPoly-1


#        if end_poly >= len(outPolyIDs):
#            end_poly = len(outPolyIDs)-1
#            print("resetting final polygon range to %d--%d" % (start_poly,end_poly))
#        if end_poly < 0:
#            end_poly = len(outPolyIDs)-1
#        nOutPolygons = end_poly - start_poly + 1
#        print("Averaging areas for %d polygons: %d to %d" %
#              (nOutPolygons, start_poly, end_poly))

        # hardwired to ensemble forcing grid formats
        times,hrus,MAP    = compAvgVal(nc_wgt,nc_in,'pcp')
        times,hrus,MATAvg = compAvgVal(nc_wgt,nc_in,'t_mean')
        times,hrus,MATRng = compAvgVal(nc_wgt,nc_in,'t_range')
        print("lengths: %d %d %d %d %d" % (len(MAP),len(MATAvg),len(MATRng),
	       len(times),len(hrus)))
        print(MAP.shape)
        writeNetCDFData(nc_out, MAP, MATAvg, MATRng, times, hrus, hru_type)

    else:
        usage()

