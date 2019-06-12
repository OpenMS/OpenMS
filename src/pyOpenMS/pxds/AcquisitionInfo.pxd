from Types cimport *
from libcpp cimport bool
from String cimport *
from MetaInfoInterface cimport *

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

