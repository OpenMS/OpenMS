from Types cimport *
from libcpp cimport bool
from MetaInfoInterface cimport *
from Acquisition cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/AcquisitionInfo.h>" namespace "OpenMS":

    cdef cppclass AcquisitionInfo(MetaInfoInterface):
        # wrap-inherits:
        #   MetaInfoInterface
        
        AcquisitionInfo() except + nogil 
        AcquisitionInfo(AcquisitionInfo &) except + nogil 

        bool operator==(AcquisitionInfo) except + nogil 
        bool operator!=(AcquisitionInfo) except + nogil 

        String getMethodOfCombination() except + nogil  # wrap-doc:Returns the method of combination
        void setMethodOfCombination(String method) except + nogil  # wrap-doc:Sets the method of combination

        Size size() except + nogil  #wrap-doc:Number a Acquisition objects
        Acquisition& operator[](size_t) except + nogil  # wrap-upper-limit:size()
        void push_back(Acquisition)  except + nogil  #wrap-doc:Append a Acquisition object
        void resize(size_t n) except + nogil 
    
        # libcpp.vector uneccessary - iteration is possible w/o, no assignment w/ or w/o
