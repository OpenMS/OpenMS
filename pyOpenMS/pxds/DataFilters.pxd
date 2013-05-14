from Types cimport *
from libcpp cimport bool
from String cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/DataFilters.h>" namespace "OpenMS":
    
    cdef cppclass DataFilters "OpenMS::DataFilters":
        DataFilters() nogil except +
        DataFilters(DataFilters) nogil except + #wrap-ignore
        Size size() nogil except +
        ### DataFilter  <](Size index) nogil except +
        ### void add(DataFilter & filter_) nogil except +
        ### void remove(Size index) nogil except +
        ### void replace(Size index, DataFilter & filter_) nogil except +
        ### void clear() nogil except +
        ### void setActive(bool is_active) nogil except +
        ### bool isActive() nogil except +
        ### bool passes(Feature & feature) nogil except +
        ### bool passes(ConsensusFeature & consensus_feature) nogil except +
        ### # TEMPLATE # bool passes(MSSpectrum[ PeakType ] & spectrum, Size peak_index) nogil except +

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/DataFilters.h>" namespace "OpenMS::DataFilters":
    cdef enum FilterType "OpenMS::DataFilters::FilterType":
        #wrap-attach:
        #    DataFilters
        INTENSITY
        QUALITY
        CHARGE
        SIZE
        META_DATA

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/DataFilters.h>" namespace "OpenMS::DataFilters":
    cdef enum FilterOperation "OpenMS::DataFilters::FilterOperation":
        #wrap-attach:
        #    DataFilters
        GREATER_EQUAL
        EQUAL
        LESS_EQUAL
        EXISTS

