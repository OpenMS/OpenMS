from Types cimport *
from libcpp cimport bool
from MetaInfoInterface cimport *
from smart_ptr cimport shared_ptr
from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string
from SpectrumSettings cimport *
from Acquisition  cimport *
from Peak1D cimport *
from String cimport *
from RangeManager cimport *
from DataArrays cimport *

cdef extern from "<OpenMS/METADATA/AcquisitionInfo.h>" namespace "OpenMS":

    cdef cppclass AcquisitionInfo(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface
        
        AcquisitionInfo()    nogil except +
        AcquisitionInfo(AcquisitionInfo)    nogil except +

        bool operator==(AcquisitionInfo) nogil except +
        bool operator!=(AcquisitionInfo) nogil except +

        String getMethodOfCombination() nogil except +
        void setMethodOfCombination(String method) nogil except +

        Size size() nogil except +
        #void reserve(size_t n) nogil except + 
        Acquisition operator[](int) nogil except + # wrap-upper-limit:size()
        #void assign(libcpp_vector[Acquisition].iterator, libcpp_vector[Acquisition].iterator) nogil except + # wrap-ignore
        libcpp_vector[Acquisition].iterator begin() nogil except +  # wrap-iter-begin:__iter__(Acquisition)
        libcpp_vector[Acquisition].iterator end()   nogil except +  # wrap-iter-end:__iter__(Acquisition)
