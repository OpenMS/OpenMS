from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from ConsensusMap cimport *
from FeatureMap cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/KERNEL/ConversionHelper.h>" namespace "OpenMS":

    cdef cppclass MapConversion:

      void convert(UInt64 input_map_index,
                   FeatureMap[Feature] input_map,
                   ConsensusMap & output_map,
                   Size n) nogil except +
    
      void convert(UInt64 input_map_index,
                   MSExperiment[Peak1D, ChromatogramPeak] & input_map,
                   ConsensusMap & output_map,
                   Size n) nogil except +
    
      void convert(ConsensusMap input_map,
                   bool keep_uids,
                   FeatureMap[Feature] & output_map) nogil except +

