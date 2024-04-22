from Types cimport *
from libcpp cimport bool
from String cimport *
from MSSpectrum cimport *
from ConsensusMap cimport *
from FeatureMap cimport *

cdef extern from "<OpenMS/PROCESSING/MISC/DataFilters.h>" namespace "OpenMS":

    cdef cppclass DataFilters "OpenMS::DataFilters":
        DataFilters() except + nogil  # compiler
        DataFilters(DataFilters &) except + nogil  # compiler 
        Size size() except + nogil 
        DataFilter operator[](size_t) except + nogil  # wrap-upper-limit:size()
        void add(DataFilter & filter_) except + nogil 
        void remove(Size index) except + nogil 
        void replace(Size index, DataFilter & filter_) except + nogil 
        void clear() except + nogil 
        void setActive(bool is_active) except + nogil 
        bool isActive() except + nogil 
        bool passes(Feature & feature) except + nogil 
        bool passes(ConsensusFeature & consensus_feature) except + nogil 
        bool passes(MSSpectrum & spectrum, Size peak_index) except + nogil 

cdef extern from "<OpenMS/PROCESSING/MISC/DataFilters.h>" namespace "OpenMS::DataFilters":

    cdef cppclass DataFilter "OpenMS::DataFilters::DataFilter":
        DataFilter() except + nogil  # compiler
        DataFilter(DataFilter &) except + nogil  # compiler
        FilterType field
        FilterOperation op
        double value
        String value_string
        String meta_name
        bool value_is_numerical
        String toString() except + nogil 
        void fromString(const String & filter_) except + nogil 
        bool operator==(DataFilter & rhs) except + nogil 
        bool operator!=(DataFilter & rhs) except + nogil 


cdef extern from "<OpenMS/PROCESSING/MISC/DataFilters.h>" namespace "OpenMS::DataFilters":
    cdef enum FilterType "OpenMS::DataFilters::FilterType":
        #wrap-attach:
        #   DataFilters
        INTENSITY
        QUALITY
        CHARGE
        SIZE
        META_DATA

cdef extern from "<OpenMS/PROCESSING/MISC/DataFilters.h>" namespace "OpenMS::DataFilters":
    cdef enum FilterOperation "OpenMS::DataFilters::FilterOperation":
        #wrap-attach:
        #   DataFilters
        GREATER_EQUAL
        EQUAL
        LESS_EQUAL
        EXISTS

