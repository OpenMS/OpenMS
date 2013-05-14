from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *

cdef extern from "<OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>" namespace "OpenMS":

    cdef cppclass PeakFileOptions:

        PeakFileOptions() nogil except +
        PeakFileOptions(PeakFileOptions) nogil except + # wrap-ignore
        
        void setMetadataOnly(bool) nogil except +
        bool getMetadataOnly()     nogil except +

        void setWriteSupplementalData(bool) nogil except +
        bool getWriteSupplementalData()     nogil except +

        void setMSLevels(libcpp_vector[int] levels) nogil except +
        void addMSLevel(Int level) nogil except +
        void clearMSLevels()       nogil except +
        bool hasMSLevels()         nogil except +
        bool containsMSLevel(int level)  nogil except +
        libcpp_vector[int] getMSLevels()    nogil except +

        void setCompression(bool) nogil except +
        bool getCompression()     nogil except +
