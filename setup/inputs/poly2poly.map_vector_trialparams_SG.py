#!/usr/bin/env python
''' Process timeseries grid into mean areal timeseries for arbitrary polygons
    Depends on mapping file from poly2poly.py and associated scripts by
    K. Sampson, NCAR '''
# Notes:  AWW: on yellowstone, needs 'module load netcdf4python'
#
#   SUMMA parameter version (for NLDAS to HUC12)
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
#  AWW Jan 2017 turned this into a parameter version (not forcing)
#               also streamlined assigning attrib from input file
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

    def getPolyWgtsIndices(self, currOutPolyID, dataIndex):
        """For the given hru id, get weight of the intersected polygons and associated
           i_index, j_index"""

        # get index corresponding to hru
        Index = [i for i, x in enumerate(self.IDmask) if x == int(currOutPolyID)]

        curr_wgt = list(self.weight[Index]) #Get wgt list for hru
        curr_intersector = list(self.intersector[Index]) #AW get intersector list for hru
        curr_intersector1 = [x-1 for x in curr_intersector] # debug SG, AW and NM Aug 01, 2019
		
        # AW now need to get matching input file 'hru' index (0-(N-1))
        # corresponding to each intersector value (which are int64 hru_id ids from
        # nldas_conus file)
        curr_input_index = list(dataIndex[curr_intersector1])

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


# ===================== main subroutine =========
def compAvgVal_scalar(wgt, var_in, varname):    
    """Compute areal weighted avg value of <varname> in <nc_in> for each hru
       based on hru's weight in <nc_wgt>"""

    print(" ")
    print("====== New Variable: %s ======= " % varname)
    print(" ")

    dataVals = var_in[:]   # extract data values from input variable object

    # ==== LOOP over output polygons
    for h in range(start_poly,(end_poly+1)):        
        # Get list of wgt, intersector index for current polygon

        (currWgtVal, currInPolyIndex) = wgt.getPolyWgtsIndices(outPolyIDs[h],dataIndex)

        nCurrCells = len(currWgtVal)  # useful limit for several loops in this block

        # ---- Calculate time series of weighted value for input polygon  ----

        # check to make sure the contributing polygons exist
        sumweights = 0
        numvoid    = 0
        validcell  = 0
        # loop through cells/polygons(ie, c) assigned to this output polygon
        #print currInPolyIndex
        print("h=%d nCurrCells=%d" % (h,nCurrCells))
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


# custom write function to handle SUMMA parameters
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
        hruId = ncfile.createVariable('hruID', 'S1', ('hru', 'strlen'))  
        hruId[:] = nc4.stringtochar(np.asarray(hrus,
                                  dtype='S{}'.format(max_strlen)))         
    else:
        # integer HRU
        hruId = ncfile.createVariable('hruID', 'i8', ('hru', ))   # edited EC
        hruId[:] = np.asarray(hrus, dtype='int')
        #hruId[:] = np.asarray(hrus, dtype='i8')        # AW for HUCs
    # add attribute    
    hruId.long_name = 'USGS HUC12 ID'
        
    # geophysical var data and attributes from input file
    print "adding data & attributes"
    ncvar = ncfile.createVariable(varname,var_out.dtype,('hru',),fill_value='-999.0')    
    ncvar[:] = var_out   # store mean areal data in netcdf file
    # for i in range(len(attNames)):
        # ncvar.setncatts({attNames[i]: attContents[i]})

    # -- from another script
    # Copy the variable attributes
    # dsout.variables[v_name].setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    # Create the variables in the file
    #for v_name, varin in dsin.variables.iteritems():
    #    dsout.createVariable(v_name, varin.datatype, varin.dimensions)
    # Copy the variables values
    #dsout.variables[v_name][:] = uuu[:]
    # -- end ... could be useful syntax --

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
               <hru_type (int or str)> <id_mapping_file.txt>
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

    if len(args) == 8:
        nc_wgts  = args[0]
        nc_infl  = args[1]
        nc_outfl = args[2]
        hru_type = args[3]
        id_map_in = args[4]  # maps from mappingfile id to datafileid
        start_poly = int(args[5])  # starts at 0 (correspond w/ netcdf conventions)
        end_poly = int(args[6])    # max possible value nOutPoly-1
        filewithvarname = args[7]     # name of input variable (just one for now)
                              # can make this a csv list to loop through
        with open(filewithvarname) as file:
            # listvars = file.readlines()
            listvars = file.read().split(',')

        for i in range(0,len(listvars)-1) :
        # for i in range(0,2) :
            varname = listvars[i]
            # SG print("The list is is %d" % (len(listvars)))
            print(varname)
                # print("checking array  .. " % (start_poly,end_poly))
        
            # read ID mapping file (translate mapping to input data index)
            # 'intersector' is index to array with data index, eg, 'hru'
            mappingId,dataIndex = np.loadtxt(id_map_in, skiprows=1,
                                             unpack=True).astype(int)
            
            # hardwired to ensemble configuration formats (hru index rather than grid)
            inPolyIDs=getInputPolyIDs(nc_infl)
            
            # get map/weight data once before starting mapping
            print("reading in weights data ... ")            
            wgts=getMappingData(nc_wgts)
            # read output polygon Ids
            outPolyIDs = wgts.getPolyIDs() # get hruID list        
            # print(outPolyIDs)
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

            # read data variable & attributes for to pass to processing and write routines
            f = nc4.Dataset(nc_infl,'r')
            var_in = f.variables[varname]
            attNames = []
            attContents = []
            attr = var_in.ncattrs()  # get attributes
            for n in range(0,len(attr)):
                attNames.extend([attr[n]])
                attContents.extend([var_in.getncattr(attr[n])])

            # maybe instead of passing all this stuff, can start w/ just writing the attrs here,
            #   then pass the data to a write function later

            # process trial params and write out in a file
            # can make var_in a list and loop over it, eventually
            hrus,var_out=compAvgVal_scalar(wgts, var_in, varname)
            # print(varname)
            # nc_outfl = nc_outfl.replace("varname",varname)
            # print(nc_outfl)
            nc_outfl_new = nc_outfl + varname + ".nc"
            writeNC_trialParam(nc_outfl_new, var_out, varname, attNames, attContents, 
                            hrus, hru_type)
            f.close() 
            
            print("The id is %d" % (i))
            
       

    else:
        usage()
