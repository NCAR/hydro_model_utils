#!/usr/bin/env python
''' Process timeseries grid into mean areal timeseries for arbitrary polygons
    Depends on mapping file from poly2poly.py and associated scripts by
    K. Sampson, NCAR '''
#
#   SUMMA attribute version (for NLDAS to HUC12)
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
#   AWW Jan 2017 * turned this into a parameter version
#                * also streamlined assigning attrib from input file
#
#                * converted to attribute version, handles both gru and hru indices
#                * also includes new 'mode' finding subroutine (for eg veg type)
#
#  AWW Mar 2017  fixed mode calculation to be area weighted
#
# =========================================================================

import sys
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

    def getPolyWgtsIndices(self, currOutPolyID, dataIndex):
        """For the given hru id, get weight of the intersected polygons and associated
           i_index, j_index"""

        # get index corresponding to hru
        Index = [i for i, x in enumerate(self.IDmask) if x == int(currOutPolyID)]

        curr_wgt = list(self.weight[Index]) #Get wgt list for hru
        curr_intersector = list(self.intersector[Index]) #AW get intersector list for hru
        curr_intersector1 = [x-1 for x in curr_intersector] # debug SG, AW and NM Aug 01, 2019
        # AW now need to get matching input forcing file 'hru' index (0-(N-1))
        # corresponding to each intersector value (which are int64 hru_id ids from
        # nldas_conus file)
        curr_input_index = list(dataIndex[(curr_intersector1)]) # debug SG, AW and NM Aug 01, 2019

        # return (self.wgt, self.i_index, self.j_index)
        return (curr_wgt, curr_input_index)

    def getPolyIDs(self):
        """ get polygon ID list of output areas (eg basins)"""
        # EC get polyid for hru list
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

def getNetCDFDim(fn, varname):
    """ Get dimension name used in <varname> from NetCDF <fn> """
    f = nc4.Dataset(fn,'r')
    dim = f.variables[varname].dimensions
    data = dict()
    for dimname in dim:
        data[dimname]=f.variables[dimname][:]
    f.close()
    return data

# AWW new subroutines ---------
def getInputPolyIDs(nc_in):
    inPolyIDs  = getNetCDFData(nc_in, 'hruId')    
    print("read input inPolyIds ('hruId') from parameter file")
    return inPolyIDs

def getMappingData(nc_wgt):
    # instantiate current wgtnc object -- creates object with mapping data
    wgts = wgtnc(nc_wgt)
    print("read weights from mapping file")
    return wgts


def get_NCvar_info(var_in):
    """ Return some variable attributes and data type"""
    attNames = []
    attContents = []
    attr = var_in.ncattrs()  # get attributes
    for n in range(0,len(attr)):
        attNames.extend([attr[n]])   # create list of names
        attContents.extend([var_in.getncattr(attr[n])])    # create list of contents
    varType = f.variables[varname].dtype

    return attNames,attContents,varType


# ===================== main subroutines ===================================

def compAvgVal(wgt, var_in, varname):    
    """Compute areal weighted avg value of <varname> in <nc_in> for each hru
       based on hru's weight in <nc_wgt>"""

    dataVals = var_in[:]   # extract data values from input variable object

    # ==== LOOP over output polygons ====
    # !! this is inefficient - weight/mapping info
    # !!   should be done once outside this call
    
    for h in range(start_poly,(end_poly+1)):        
        # Get list of wgt, intersector index for current polygon

        (currWgtVal, currInPolyIndex) = wgt.getPolyWgtsIndices(outPolyIDs[h],dataIndex)

        nCurrCells = len(currWgtVal)  # useful limit for several loops in this block

        print currInPolyIndex
        print currWgtVal
        print("h=%d nCurrCells=%d" % (h,nCurrCells))

        # check to make sure the contributing polygons exist
        # loop through cells/polygons (ie, c below) assigned to this output polygon
        sumweights = 0
        numvoid    = 0
        validcell  = 0
        for c in range(0, nCurrCells):
            validcell = validcell + 1

            # check if inPoly at index is missing or it's present but the first
            #   value is no_data -- if either then set weight to zero
            if currInPolyIndex[c] < 0 or dataVals[currInPolyIndex[c]] > 1.e+19:
                currWgtVal[c] = 0
                numvoid   = numvoid + 1

            sumweights = sumweights + currWgtVal[c]

        print(" var %s outPoly %d %s: num voids = %d out of %d" %
              (varname, h, outPolyIDs[h], numvoid, validcell))
        print(" new sum of weights is %f before normalizing" % sumweights)

        # if there is contributing data to average and we can proceed...
        if sumweights > 0:
            # renormalize weights
            for c in range(0, nCurrCells):
                currWgtVal[c] = (currWgtVal[c]/sumweights)

            oriVals = 0  # create zero valu to start array for sums
            
            # now apply each wgt to cell data for all times
            for c in range(0, nCurrCells):

                # if weight is above zero, add data from new cell
                if currWgtVal[c]>0:
                    # apply weight values
                    Vals1 = currWgtVal[c]*dataVals[currInPolyIndex[c]]
                    oriVals = np.append(oriVals, Vals1)

            # now sum up weighted input polygon cells, delete temporary array
            Vals2=oriVals.sum()

            # now store time vector of weighted average results
	    if 'wgtedVals' in locals():
                wgtedVals=np.append(wgtedVals, Vals2)
                activeOutPolyIDs = np.append(activeOutPolyIDs, outPolyIDs[h])
            else:
  	        # first time, create the start fof the list to be summed
                wgtedVals = Vals2
                activeOutPolyIDs = outPolyIDs[h]                

            
        # else none of the contributing cells are present and avg can't be formed
        else:
            print("WARNING: no valid input contributing areas found")
            print("WARNING: all cells (%d of %d)" % (numvoid, validcell) +
                  " in hru %d are void or missing: omitting HRU from output" % (h))

        print("--------------------")

    return activeOutPolyIDs,wgtedVals


def findModeVal(wgt, var_in, varname):    
    """Compute mode value of <varname> in <nc_in> for each hru based on hru's
       weight in <nc_wgt>"""

    dataVals = var_in[:]   # extract data values from input variable object

    # ==== LOOP over output polygons
    for h in range(start_poly,(end_poly+1)):        
        # Get list of wgt, intersector index for current polygon

        (currWgtVal, currInPolyIndex) = wgt.getPolyWgtsIndices(outPolyIDs[h],dataIndex)

        nCurrCells = len(currWgtVal)  # useful limit for several loops in this block

        print currInPolyIndex
        print("h=%d nCurrCells=%d" % (h,nCurrCells))

	# set up counter list to calculate weighted mode (start w/ 0s)
	typeSum = [0]*100   # use 100 as max possible mode type index (eg vegtype)

        for c in range(0, nCurrCells):
            if currInPolyIndex[c] > 0:
                #print("%d %d" % (c, dataVals[currInPolyIndex[c]]))
		if dataVals[currInPolyIndex[c]] > (len(typeSum)-1):
			print("type index exceeds range for counter in findModeVal()",
				"increase range in code")
			sys.exit(2)
		else:
			typeSum[dataVals[currInPolyIndex[c]]] += currWgtVal[c]

	if(max(typeSum)) > 0:
	    mode = [i for i, x in enumerate(typeSum) if x == max(typeSum)]
            print(" found %s mode data %f" % (varname, mode[0]))  # take 1st if multi modes
	
            # store weighted average results
	    if 'wgtedVals' in locals():
                wgtedVals=np.append(wgtedVals, mode)
                activeOutPolyIDs = np.append(activeOutPolyIDs, outPolyIDs[h])
            else:
  	        # first time, create the start fof the list to be summed
                wgtedVals = mode
                activeOutPolyIDs = outPolyIDs[h]                

            
        # else none of the contributing cells are present and avg can't be formed
        else:
            print("WARNING: no valid input contributing data found")
            print("WARNING: omitting hru = %d from output" % (h))

        print("--------------------")

    return activeOutPolyIDs,wgtedVals


def findActiveHRUs(wgt, var_in, varname):    
    """Compute areal weighted avg value of <varname> in <nc_in> for each hru
       based on hru's weight in <nc_wgt>"""

    dataVals = var_in[:]   # extract data values from input variable object

    # ==== LOOP over output polygons ====
    # !! this is inefficient - weight/mapping info
    # !!   should be done once outside this call
    
    for h in range(start_poly,(end_poly+1)):        
        # Get list of wgt, intersector index for current polygon
        #   don't need weights for this actually

        (currWgtVal, currInPolyIndex) = wgt.getPolyWgtsIndices(outPolyIDs[h],dataIndex)

        nCurrCells = len(currWgtVal)  # useful limit for several loops in this block

        print currInPolyIndex
        print("h=%d nCurrCells=%d" % (h,nCurrCells))

        # check to make sure the contributing polygons exist
        # loop through cells/polygons (ie, c below) assigned to this output polygon
        sumweights = 0
        numvoid    = 0
        validcell  = 0
        for c in range(0, nCurrCells):
            validcell = validcell + 1

            # check if inPoly at index is missing or it's present but the first
            #   value is no_data -- if either then set weight to zero
            if currInPolyIndex[c] < 0 or dataVals[currInPolyIndex[c]] > 1.e+19:
                currWgtVal[c] = 0
                numvoid   = numvoid + 1

            sumweights = sumweights + currWgtVal[c]

        print(" var %s outPoly %d %s: num voids = %d out of %d" %
              (varname, h, outPolyIDs[h], numvoid, validcell))
        print(" new sum of weights is %f before normalizing" % sumweights)

        # if there is contributing data to average and we can proceed...
        if sumweights > 0:

            # now store time vector of weighted average results
	    if 'activeOutPolyIDs' in locals():
                activeOutPolyIDs = np.append(activeOutPolyIDs, outPolyIDs[h])
            else:
  	        # first time, create the start fof the list to be summed
                activeOutPolyIDs = outPolyIDs[h]                
            
        # else none of the contributing cells are present and avg can't be formed
        else:
            print("WARNING: no valid input contributing areas found")
            print("WARNING: all cells (%d of %d)" % (numvoid, validcell) +
                  " in hru %d are void or missing: omitting HRU from output" % (h))

        print("--------------------")

    print("number of activeOutPolyIDs is %d" % len(activeOutPolyIDs))

    return activeOutPolyIDs


# ======= writing subroutines ===========

# write dimensions and dimension variables to netcdf output file
def writeNC_dims(fn,  hrus, hru_type):    
    """ Write <vars>[hru] array in netCDF4 file,<fn> and variable of
        <varname> """

    print "starting output file"
    nc_out = nc4.Dataset(fn, 'w', format='NETCDF4')

    # Create dimensions
    dim_hru = nc_out.createDimension('hru', len(hrus))
#    dim_gru = nc_out.createDimension('gru', len(hrus))    

    # --- Create HRU and GRU ID variables (can be either int or string)
    # --- currently these are all the same
    if hru_type == 'str':
        # string HRU (need to add string length)
        max_strlen = 20 
        dim_str = nc_out.createDimension('strlen', max_strlen)
        hruId = nc_out.createVariable('hruId', 'S1', ('hru', 'strlen'))  
        hruId[:] = nc4.stringtochar(np.asarray(hrus,
                                  dtype='S{}'.format(max_strlen)))
#        gruId = nc_out.createVariable('gruId', 'S1', ('hru', 'strlen'))  
#        gruId[:] = nc4.stringtochar(np.asarray(hrus,
#                                  dtype='S{}'.format(max_strlen)))
        hru2gruId = nc_out.createVariable('hru2gruId', 'S1', ('hru', 'strlen'))  
        hru2gruId[:] = nc4.stringtochar(np.asarray(hrus,
                                  dtype='S{}'.format(max_strlen)))        

    else:
        # integer HRUs
        hruId = nc_out.createVariable('hruId', 'i8', ('hru', ))  
        hruId[:] = np.asarray(hrus, dtype='int')
#        gruId = nc_out.createVariable('gruId', 'i8', ('gru', ))  
#        gruId[:] = np.asarray(hrus, dtype='int')
        hru2gruId = nc_out.createVariable('hru2gruId', 'i8', ('hru', ))  
        hru2gruId[:] = np.asarray(hrus, dtype='int')        

    # add attributes
    hruId.long_name = 'Ids for hydrological response units'
#    gruId.long_name = 'Id of group of response unit (GRU) -- USGS HUC12 ID'
    hru2gruId.long_name = 'Id of GRU to which the HRU belongs'    

    return nc_out
    # leave netcdf file open


# append variables to netcdf output file
def append_NC_vars(nc_out, newVarName, newVarDim, newVarType, newVarVals,
                   attNames, attContents):

    """ Write <vars>[hru] array in netCDF4 file,<fn> and variable of
        <varname> """

    ncvar = nc_out.createVariable(newVarName, newVarType, (newVarDim,))    
    ncvar[:] = newVarVals   # store data in netcdf file

    # add attributes that were in input file
    for i in range(len(attNames)):
        if(attNames[i] != '_FillValue'):   # gives problems with one variable
            ncvar.setncatts({attNames[i]: attContents[i]})
            print("setting %s: %s" % (attNames[i], attContents[i]))
            

# custom write function to handle SUMMA forcing variables
def writeNC_trialParam(fn,  var_out, varname, attNames, attContents, hrus, hru_type):    
    """ Write <vars>[hru] array in netCDF4 file,<fn> and variable of
        <varname> """

    print "writing output file"
    ncfile = nc4.Dataset(fn, 'w', format='NETCDF4')

    # Data should be 1D [hru]
    dimSize = var_out.shape[0]
    dim_1 = ncfile.createDimension('hru', dimSize)  # hru axis

    # ==== Populate netcdf file variables =====
    
    # --- Create HRU ID variable (can be either int or string)
    if hru_type == 'str':
        # string HRU (need to add string length)
        max_strlen = 20  # EC
        dim_str = ncfile.createDimension('strlen', max_strlen)
        hruId = ncfile.createVariable('hru', 'S1', ('hru', 'strlen'))  
        hruId[:] = nc4.stringtochar(np.asarray(hrus,
                                  dtype='S{}'.format(max_strlen)))         
    else:
        # integer HRU
        hruId = ncfile.createVariable('hru', 'i8', ('hru', ))   # edited EC
        hruId[:] = np.asarray(hrus, dtype='int')
        #hruId[:] = np.asarray(hrus, dtype='i8')        # AW for HUCs
    # add attribute    
    hruId.long_name = 'USGS HUC12 ID'
        
    # geophysical var data and attributes from input file
    print "adding data & attributes"
    ncvar = ncfile.createVariable(varname,var_out.dtype,('hru',),fill_value='-9999.0')    
    ncvar[:] = var_out   # store mean areal data in netcdf file
    for i in range(len(attNames)):
        ncvar.setncatts({attNames[i]: attContents[i]})

    # Write basic global attribute
    ncfile.history = 'Created ' + time.ctime(time.time())
    ncfile.source = os.path.dirname(os.path.abspath(__file__))+__file__[1:]

    ncfile.close()
    # end writing netcdf output function

    


############################################
#                  Main                    #
############################################
use = '''
Usage: %s -[h] <weight_netCDF> <input_netCDF> <output_netCDF>
               <hru_type (int or str)> <id_mapping_file.txt> <start_poly> <end_poly>
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
        nc_wgts_fname  = args[0]
        nc_in_fname  = args[1]   # attributes file to be mapped
        nc_out_fname = args[2]   # new attributes file
        hru_type = args[3]   # int or str
        id_map_fname = args[4]  # maps from mappingfile id to datafileid
        start_poly = int(args[5])  # starts at 0 (correspond w/ netcdf conventions)
        end_poly = int(args[6])    # max possible value nOutPoly-1

        # read ID mapping file (translate mapping to input data index)
        # 'intersector' is index to array with data index, eg, 'hru'
        mappingId,dataIndex = np.loadtxt(id_map_fname, skiprows=1,
                                         unpack=True).astype(int)

        # dataIndex = np.where(dataIndex==-99, -99, dataIndex-1)
        
        # hardwired to ensemble forcing formats (hru index rather than grid)
        inPolyIDs=getInputPolyIDs(nc_in_fname)

        # get map/weight data once before starting mapping
        # right now  do this for all even though only writing some
        #   can make more efficient by subsetting (doesn't take long though)
        print("reading in weights data ... ")            
        wgts=getMappingData(nc_wgts_fname)
        # read output polygon Ids
        outPolyIDs = wgts.getPolyIDs() # get hruID list        
        (wgtsAll,outIDmask, intersectorsAll) = wgts.getWgtDataAll()
        print(" ")

        # calculate/reset run limits for current script call, dep. on data
        nOutPolygons = len(outPolyIDs)

        if end_poly >= len(outPolyIDs):
            end_poly = len(outPolyIDs)-1
            print("resetting final polygon range to %d--%d" % (start_poly,end_poly))
        if end_poly < 0:
            end_poly = len(outPolyIDs)-1
        nOutPolygons = end_poly - start_poly + 1
	print("Averaging areas for %d polygons: %d to %d" %
              (nOutPolygons, start_poly, end_poly))

        # set input file
        f = nc4.Dataset(nc_in_fname,'r')

        # now work on variables that will be averaged (all indexed by hru)
        #   could read these from a dictionary to determine how to treat them
        for varname in ("contourLength", "tan_slope", "mHeight", "HRUarea", "elevation",
                        "latitude", "longitude"):

#	    break

            # when processing the first variable, get hrus and start output file

            # -- needs update because this can write contourLength rather than assess it twice
            if varname == "contourLength":
                # find active HRUs (those that will have valid data
                var_in = f.variables["contourLength"]        
                hrus = findActiveHRUs(wgts, var_in, "contourLength")

                # initialize netcdf file by storing dimensions and hru variables
                nc_out = writeNC_dims(nc_out_fname, hrus, hru_type)                
        
            # read data variable & attributes for to pass to processing and write routines 
            print("\n====== New Variable: %s =======" % varname)           
            var_in = f.variables[varname]

            # areal average attributes
            hrus,var_out = compAvgVal(wgts, var_in, varname)

            # get variable attributes
            attNames,attContents,varType = get_NCvar_info(var_in)

            # add to output netcdf file
            print("adding variable %s to output file" % varname)
            append_NC_vars(nc_out, varname, 'hru', varType, var_out, attNames, attContents)
            
        # special case variable that is all zeros for now
        for varname in ("downHRUindex",):

#	    break

            # read data variable & attributes for to pass to processing and write routines
            print("\n====== New Variable: %s =======" % varname)
            var_in = f.variables[varname]

            # set variable to zero, reuse 'hrus' from earlier calls
            var_out = np.zeros((1,len(hrus)))

            # get variable attributes
            attNames,attContents,varType = get_NCvar_info(var_in)

            # add to output netcdf file
            print("adding variable %s to output file" % varname)
            append_NC_vars(nc_out, varname, 'hru', varType, var_out, attNames, attContents)

        # variables that aren't averaged but take mode
        for varname in ("slopeTypeIndex","soilTypeIndex","vegTypeIndex"):

# testing
#        for varname in ("vegTypeIndex",):

  	    # find active HRUs (those that will have valid data)
 #           var_in = f.variables["contourLength"]
 #           hrus = findActiveHRUs(wgts, var_in, "contourLength")

            # initialize netcdf file by storing dimensions and hru variables
 #           nc_out = writeNC_dims(nc_out_fname, hrus, hru_type)


 
            # read data variable & attributes for to pass to processing and write routines
            print("\n====== New Variable: %s =======" % varname)
            var_in = f.variables[varname]

            # set variable to zero, reuse 'hrus' from earlier calls
            hrus,var_out = findModeVal(wgts, var_in, varname)            

            # get variable attributes
            attNames,attContents,varType = get_NCvar_info(var_in)

            # add to output netcdf file
            print("adding variable %s to output file" % varname)
            append_NC_vars(nc_out, varname, 'hru', varType, var_out, attNames, attContents)
        
        
        nc_out.close()        
        f.close()        

    else:
        usage()


