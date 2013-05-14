from Types cimport *
from libcpp cimport bool
from String cimport *
from MSSpectrum cimport *
from DataFilters cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/DataFilters.h>" namespace "OpenMS::DataFilters":
    
    cdef cppclass DataFilter "OpenMS::DataFilters::DataFilter":
        DataFilter() nogil except +
        DataFilter(DataFilter) nogil except + #wrap-ignore
        FilterType field
        FilterOperation op
        DoubleReal value
        String value_string
        String meta_name
        bool value_is_numerical
        String toString() nogil except +
        void fromString(String & filter_) nogil except +
        bool operator==(DataFilter & rhs) nogil except +
        bool operator!=(DataFilter & rhs) nogil except +

