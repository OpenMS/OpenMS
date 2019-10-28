from Types cimport *
from libcpp cimport bool
from String cimport *
from MetaInfoInterface cimport *
from Acquisition cimport *

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

        Size size() nogil except + #wrap-doc:Number a Acquisition objects
        Acquisition operator[](int) nogil except + # wrap-upper-limit:size()
        void push_back(Acquisition)  nogil except + #wrap-doc:Append a Acquisition object
        void resize(size_t n) nogil except +

