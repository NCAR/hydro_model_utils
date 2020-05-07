#!/usr/bin/env python
''' Process timeseries grid into mean areal timeseries for arbitrary polygons
    Depends on mapping file from poly2poly.py and associated scripts by
    K. Sampson, NCAR '''
# Notes:  AWW: on yellowstone, needs 'module load netcdf4python'
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
#
#   AWW Jan 2017 substantial rewrite:  moved some repetitive data reads out of
#                  loops, cleaned up variable names to add clarity
#                  added an ID-mapping capability in case the data array indices
#                  doesn't match the hru-id mapping (as in subsetted domains)
#                  configured for vector rather than grid mapping, and
#                  added subsetting of hru-domain for testing/development
#                  variables match the full NLDAS hourly forcings (7 vars)
#
# =========================================================================

import sys
#sys.path.append("/d6/nijssen/anaconda/lib/python2.7/site-packages/")

import os
import time
import getopt
import numpy as np
import netCDF4 as nc4
#import xarray as xr

############################################
#              Class                       #
############################################
class wgtnc:
    # """object of basin netCDF including hru and areal weight/i_index/j_index
    """AW object of basin netCDF including hru and areal weight/intersector
           of others polygons  """
    def __init__(self,nc):
        """Initialization """
        self.nc=nc

    def getWgtDataAll(self):
        """AW: Read the weighting file once, return the key data   """
        self.weight      = getNetCDFData(self.nc, 'weight')
        self.IDmask      = getNetCDFData(self.nc, 'IDmask')
        self.intersector = getNetCDFData(self.nc, 'intersector')
        print("read mapping file data")
        
        return self.weight,self.IDmask,self.intersector


    def matchIntersectorToInputIndex(self, dataIndex):
        """AW: match up the IDs to the input variable indices
           ideally they would match, but in some cases they don't
           this is now combined with the next subroutine -- not used"""

        # never got this to work fast enough (using enumerate), so now read in a
        #   pre-prepared id map file instead

        # duplicate variable, then repopulate
        self.input_index = self.intersector
        self.input_index = dataIndex[self.intersector]        
        print("matched mapping file intersectors to data indexes")
        
        return self.input_index


    def getPolyWgtsIndices(self, currOutPolyID):
        """For the given hru id, get weight of the intersected polygons and associated
           i_index, j_index"""

        # get indices (record #s in vector) corresponding to hru matches
        # note that currOutPolyID is a string ... IDmask is an integer
 
        Index = [i for i, x in enumerate(self.IDmask) if x == int(currOutPolyID)]
        print("currpoly = %s " % currOutPolyID )
        #print Index

        # use the indices to pull out the matching weights and hruID intersection
        curr_wgt = list(self.weight[Index]) #Get wgt list for hru
        curr_intersector = list(self.intersector[Index]) #AW get intersector list for hru
        curr_intersector1 = [x-1 for x in curr_intersector] # debug SG, AW and NM Aug 01, 2019
        # AW now find the data vector indices corresponding to each intersector value 
        # (which are int64 hru_id ids from nldas_conus file)
        curr_input_index = list(dataIndex[curr_intersector1])

        return (curr_wgt, curr_input_index)

    

    def getWgtHru(self, currOutPolyID, inPolyIDs ):
        """For the given hru id, get weight of the intersected polygons and associated
           dataIndex"""

        # AWW NOT USED NOW (replaced by above)

        # !AWW: does this really READ the mapping netcdf with every single HRU call?
        
        wgtAll = getNetCDFData(self.nc, 'weight')
        polyAll = getNetCDFData(self.nc, 'IDmask')
        # EC get i_index --- LON direction
        #i_indexAll = getNetCDFData(self.nc, 'i_index')
        # EC get j_index --- LAT direction
        #j_indexAll = getNetCDFData(self.nc, 'j_index')
        # AW get intersector value (ID of starting polygon/'grid')
        #   used in this version
        intersectorAll = getNetCDFData(self.nc, 'intersector')

        #PolyList = self.getPolyIDs() # EC-get hru id list (used?
        
        # EC-get index corresponding to hru
        Index = [i for i, x in enumerate(polyAll) if x == currOutPolyID]
        print("Index:")
        print Index
        sys.exit()  # should not be here anyway

        
        self.wgt = list(wgtAll[Index]) #Get wgt list for hru
        #self.i_index = list(i_indexAll[Index]) #Get i_index list for hru
        #self.j_index = list(j_indexAll[Index]) #Get j_index list for hru
        self.intersector = list(intersectorAll[Index]) #AW get intersector list for hru

        # AW now need to get matching input forcing file 'hru' index (0-(N-1))
        # corresponding to each intersector value (which are int64 hru_id ids from
        # nldas_conus file)
        #self.intersector_index = self.intersector
        #print("%s %d\n" % (self.intersector[1], self.intersector_index[1]))

        #is this a good place to swap in the incremental 'hru'?  is that needed?

        # return (self.wgt, self.i_index, self.j_index)
        return (self.wgt, self.intersector)  # AW make this the index, not intersector


    def getPolyIDs(self):
        """ get polygon ID list of output areas (eg basins)"""
        # EC get polyid for hru list
        # note these are strings
        self.polyid = list(getNetCDFData(self.nc, 'polyid'))
        return self.polyid

# --- end of class wgtnc definition -----

########################################################################
#                                Modules                               #
########################################################################

def getNetCDFData(fn, varname):
    """Read <varname> variables available to be mapped from NetCDF <fn> """
    f = nc4.Dataset(fn,'r')
    data = f.variables[varname][:]
    f.close()

#    ds = xr.open_dataset(fn)
#    data = ds[varname]

    return data
    
# next fcn not used
def getNetCDFVarDimNames(fn, varname):
    """ Get dimension names used in <varname> from NetCDF <fn> """
    f = nc4.Dataset(fn,'r')
    dim = f.variables[varname].dimensions

    data = dict()
    for dimname in dim:
        data[dimname]=f.variables[dimname][:]
    f.close()
    return data

# custom write function to handle SUMMA forcing variables
def writeNetCDFData(fn,  prec, temp, pres, dslr, lwav, shum, wind,
                    times, hrus, hru_type):    
    """ Write <vars>[time,hru] array in netCDF4 file,<fn> and variable of
        <varname> """
    ncfile = nc4.Dataset(fn, 'w', format='NETCDF4')

    # Data should be 2D [time x hru]
    dim2size = prec.shape[1]
    # time axis - record dimension
    dim_1 = ncfile.createDimension('time', None)
    dim_2 = ncfile.createDimension('hru', dim2size)  # hru axis

    # ==== Populate netcdf file variables =====
    # Create netcdf file Variables and define dimensions variable to hold the var
    
    # --- Create HRU ID variable (can be either int or string)
    if hru_type == 'str':
        # string HRU
        max_strlen = 20  # EC
        # EC: added string length for hrus
        dim_3 = ncfile.createDimension('strlen', max_strlen)
        hruId = ncfile.createVariable('hruId', 'S1', ('hru', 'strlen'))  
        hruId[:] = nc4.stringtochar(np.asarray(hrus,
                                  dtype='S{}'.format(max_strlen)))         
    else:
        # integer HRU
        hruId = ncfile.createVariable('hruId', 'i8', ('hru', ))   # edited EC
        #hruId[:] = np.asarray(hrus, dtype='int')
        hruId[:] = np.asarray(hrus, dtype='i8')        # AW for HUCs
    
    #hruId.long_name = 'USGS HUC12 ID'
    hruId.long_name = 'USGS ID'  # for camels
        
        
    # --- Times ---
    # t = ncfile.createVariable('time', 'i4', ('time', )) # int for seconds
    t = ncfile.createVariable('time', 'f8', ('time', ))   # float for 'days' # SG July 17, 2019, converted to f8 from f4
    t[:] = times   # add time data    
    #t.time_origin = '1990-JAN-01 00:00:00'
    #t.title = 'time'
    t.long_name = 'Observation time'
    t.units = 'days since 1990-01-01 00:00:00'        
    t.calendar = 'standard'
    #t.axis = 'T'    

    # geophysical vars: pptrate, LWRadAtm, SWRadAtm, airpres, airtemp, spechum, windspd
    pptr = ncfile.createVariable('pptrate','f4',('time','hru',),fill_value='-999.0')
    pptr[:,:] = prec   # for mean area precip, etc.
    pptr.long_name = 'Precipitation rate'
    pptr.units = 'kg m-2 s-1'
    pptr.setMissing = '1.e+20'
    #pptr.associate = 'time hru'    
    #map.axis = 'TX'

    airt = ncfile.createVariable('airtemp','f4',('time','hru',),fill_value='-999.0')
    airt[:,:] = temp  
    airt.long_name = 'air temperature at the measurement height'
    airt.units = 'K'
    airt.setMissing = '1.e+20'

    airp = ncfile.createVariable('airpres','f4',('time','hru',),fill_value='-999.0')
    airp[:,:] = pres  
    airp.long_name = 'air pressure at the measurement height'
    airp.units = 'Pa'
    airp.setMissing = '1.e+20'

    sphu = ncfile.createVariable('spechum','f4',('time','hru',),fill_value='-999.0')
    sphu[:,:] = shum  
    sphu.long_name = 'specific humidity at the measurement height'
    sphu.units = 'g g-1'
    sphu.setMissing = '1.e+20'
    
    wspd = ncfile.createVariable('windspd','f4',('time','hru',),fill_value='-999.0')
    wspd[:,:] = wind
    wspd.long_name = 'wind speed at the measurement height'
    wspd.units = 'm s-1'
    wspd.setMissing = '1.e+20'

    lrad = ncfile.createVariable('LWRadAtm','f4',('time','hru',),fill_value='-999.0')
    lrad[:,:] = lwav  
    lrad.long_name = 'downward longwave radiation at the upper boundary'
    lrad.units = 'W m-2'
    lrad.setMissing = '1.e+20'
    
    srad = ncfile.createVariable('SWRadAtm','f4',('time','hru',),fill_value='-999.0')
    srad[:,:] = dslr
    srad.long_name = 'downward shortwave radiation at the upper boundary'
    srad.units = 'W m-2'
    srad.setMissing = '1.e+20'    
    
    # martyn's data_step variable for summa
    tstep = ncfile.createVariable('data_step','d',dimensions=())
    tstep[:] = 3600.0   # assign value; should get from input
    # tstep[:] = 10800.0   # assign value; should get from input
    tstep.long_name = 'data step length in seconds'
    tstep.units = 'seconds' 

    # Write basic global attribute
    ncfile.history = 'Created ' + time.ctime(time.time())
    ncfile.source = os.path.dirname(os.path.abspath(__file__))+__file__[1:]

    ncfile.close()

    # end writing netcdf output function

# AWW new subroutines ---------
def getInputDataDims(nc_in):
    timeVals   = getNetCDFData(nc_in,'time')
    inPolyIDs  = getNetCDFData(nc_in, 'hruId')    
    #timeDim   = timeVals.shape[0]      # time dimension
    #inPolyDim = InputHRUsAll.shape[0]  # input data dimension
    print("read time and inPolyIds ('hruId') from forcing file")

    return timeVals,inPolyIDs


def getMappingData(nc_wgt):
    # instantiate current wgtnc object -- creates object with mapping data
    wgts = wgtnc(nc_wgt)
    print("read weights from mapping file")
    
    return wgts


# ==================================== main subroutine =========
#def compAvgVal(nc_wgt, nc_in, varname):
def compAvgVal(wgt, nc_in, varname):
    """Compute areal weighted avg value of <varname> in <nc_in> for each hru
       based on hru's weight in <nc_wgt>"""

    print(" ")
    print("====== New Variable: %s ======= " % varname)
    print(" ")

    # instantiate current wgtnc object -- creates object with mapping data
    #wgt = wgtnc(nc_wgt)  # moved out so it doesn't repeat
    outPolyIDs = wgt.getPolyIDs() # get hruID list (these are strings)
    #print("polyIDs[0] = %s" % polyIDs[0])    
    
    # now read the input timeseries file
    print("reading input timeseries data file")

    dataVals  = getNetCDFData(nc_in, varname)     # Get data values
    # could speed this up by not importing time after the first time
    timeVals  = getNetCDFData(nc_in,'time')
    inPolyIDs = getNetCDFData(nc_in, 'hruId')    

#    print dataVals.shape
#    print dataVals.values[:,1]
#    sys.exit()

    timeDim   = timeVals.shape[0]      # time dimension
    inPolyDim = dataVals.shape[1]  # input data dimension

    # read in all the weight data before loop through IDs
    (wgtsAll,outIDmask, intersectorsAll) = wgt.getWgtDataAll()
    #print("intersectorsAll[0] = %d" % intersectorsAll[0])
    print(" ")

    # ==== LOOP over output polygons
    for h in range(start_poly,(end_poly+1)):        
        # Get list of wgt, intersector index for current polygon

        (currWgtVal, currInPolyIndex) = wgt.getPolyWgtsIndices(outPolyIDs[h])

        nCurrCells = len(currWgtVal)  # useful limit for several loops in this block

        # ---- Calculate time series of weighted value for input polygon  ----

        # check to make sure the contributing polygons exist
        sumweights = 0
        numvoid    = 0
        validcell  = 0
        # loop through cells/polygons(ie, c) assigned to this output polygon
        for c in range(0, nCurrCells):
            validcell = validcell + 1

            # check if inPoly at index is missing or it's present but the first
            #   value is no_data -- if either then set weight to zero
            if currInPolyIndex[c] < 0 or dataVals[0,currInPolyIndex[c]] > 1.e+19:
                currWgtVal[c] = 0
                numvoid   = numvoid + 1

            sumweights = sumweights + currWgtVal[c]

        print(" var %s outPoly %d %s: num voids = %d out of %d" % (varname, h,
                                                               outPolyIDs[h],
                                                               numvoid,
                                                               validcell))
        print(" new sum of weights is %f before normalizing" % sumweights)

        # if there is contributing data to average and we can proceed...
        if sumweights > 0:
            # renormalize weights
            for c in range(0, nCurrCells):
                currWgtVal[c] = (currWgtVal[c]/sumweights)

            # create zero vector to start array for sums
            oriVals = 0 * dataVals[:,0].reshape(-1,1)
            
            # now apply each wgt to cell data for all times
            for c in range(0, nCurrCells):

                # if weight is above zero, add data from new cell
                if currWgtVal[c]>0:
                    # apply weight values -- this is a vector with time
                    Vals1 = currWgtVal[c]*dataVals[:,currInPolyIndex[c]].reshape(-1,1)
                    #Val1=wgtval[c]*dataVal[:,yj,xi].reshape(-1,1)
                    oriVals = np.concatenate((oriVals, Vals1), axis=1)


            # now sum up weighted input polygon cells, delete temporary array
            Vals2=oriVals.sum(axis=1).reshape(-1,1)
            #print Vals2.shape
            #print oriVals.shape
            del oriVals

            # now store time vector of weighted average results
	    if 'wgtedVals' in locals():
                wgtedVals=np.concatenate((wgtedVals, Vals2), axis=1)
                activeOutPolyIDs = np.append(activeOutPolyIDs, outPolyIDs[h])            
            else:
  	        # first time, create the start fof the multi-column array
                wgtedVals = Vals2
                activeOutPolyIDs = outPolyIDs[h]                


        # else none of the contributing cells are present and avg can't be formed
        else:
            print("WARNING: no valid input contributing areas found")
            print("WARNING: all cells (%d of %d)" % (numvoid, validcell) +
                  " in hru %d are void or missing: omitting HRU from output" % (h))

		#  " in hru %d are void or missing: outputting voids" % (h))

            # could also assign voids and move the concatenate above down below
            #for tm in range(len(Vals2)):
                #Vals2[tm] = 1.e+20  # orig
		#Vals2 = (0 * dataVals[:,0] + 1.e+20).reshape(-1,1)                


        print("--------------------")

# NOTE: this could fail if no valid weighted avgs were made, and wgtedVals was not created

    return timeVals,activeOutPolyIDs,wgtedVals
    #return timeVals,outPolyIDs,wgtedVals

############################################
#                  Main                    #
############################################
use = '''
Usage: %s -[h] <weight_netCDF> <input_netCDF> <output_netCDF>
               <hru_type (int or str)> <id_mapping_file.txt>
               <start_out_rec=0> <end_out_rec>
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
            raise OptionError(opt)
            usage()

    if len(args) == 7:
        nc_wgts    = args[0]
        nc_in      = args[1]
        nc_out     = args[2]
        hru_type   = args[3]
        id_map_in  = args[4]
        start_poly = int(args[5])  # starts at 0 (correspond w/ netcdf conventions)
        end_poly   = int(args[6])    # max possible value nOutPoly-1        


        # read ID mapping file (translate mapping to input data index)
        a,b = np.loadtxt(id_map_in, skiprows=1, unpack=True).astype(int)        
        mappingId = a.astype(int)
        dataIndex = b.astype(int)
        
        # hardwired to ensemble forcing formats (hru index rather than grid)
        times,inPolyIDs=getInputDataDims(nc_in)

        # get map/weight data once before starting mapping
        wgts=getMappingData(nc_wgts)
        # read output polygon Ids
        outPolyIDs = wgts.getPolyIDs() # get hruID list        
        print(" ")        

        # calculate run limits
        nOutPolygons = len(outPolyIDs)
        if end_poly >= len(outPolyIDs):
            end_poly = len(outPolyIDs)-1
            print("resetting final polygon range to %d--%d" % (start_poly,end_poly))
        if end_poly < 0:
            end_poly = len(outPolyIDs)-1
        nOutPolygons = end_poly - start_poly + 1
	print("Averaging areas for %d polygons: %d to %d" %
              (nOutPolygons, start_poly, end_poly))

        # compute average values for first variable
        times,hrus,prec=compAvgVal(wgts, nc_in, 'pptrate')
        times,hrus,temp=compAvgVal(wgts, nc_in, 'airtemp')
        times,hrus,pres=compAvgVal(wgts, nc_in, 'airpres')        
        times,hrus,dslr=compAvgVal(wgts, nc_in, 'SWRadAtm')
        times,hrus,lwav=compAvgVal(wgts, nc_in, 'LWRadAtm')
        times,hrus,shum=compAvgVal(wgts, nc_in, 'spechum')
        times,hrus,wind=compAvgVal(wgts, nc_in, 'windspd')                                
        writeNetCDFData(nc_out, prec, temp, pres, dslr, lwav, shum, wind,
                        times, hrus, hru_type)        

    else:
        usage()
