# S. Frank, Realta Fusion, 2023/08
#
# Helper file containing cql3d output file reader 

import numpy as np
import scipy.io as spio

class cql3d_ncdf():
    def __init__(self,cql_mnemonic='cql3d',cql_tag=None):

        #open cql3d netcdf
        try:
            if cql_tag == None:
                self.cql_tagname = cql_mnemonic
            else:
                self.cql_tagname = cql_mnemonic+cql_tag
                
            self.cql_filename = self.cql_tagname+".nc"
                
            cql_nc = spio.netcdf_file(self.cql_filename,'r')
        except:
            print('cql3d_file initialization failed: '\
                  'could not find ncdf: ',self.cql_filename)
            raise Exception('cql3d_file initialization failed: could '\
                            'not find ncdf: ',self.cql_filename)

        #read in cdf dimension
        #creates an empty dict. and fills it up with the cdf dims.
        self.dim = {}
        for name in cql_nc.dimensions.keys():
            self.dim[name] = \
                np.copy(cql_nc.dimensions[name])

        #read in the cdf variables
        #creates an empty dict. and fills it up with the cdf vars.
        self.var = {}
        for name in cql_nc.variables.keys():
            self.var[name] = \
                np.copy(cql_nc.variables[name].data)

        #close the cdf file
        cql_nc.close()
