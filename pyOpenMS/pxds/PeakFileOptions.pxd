from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from DRange cimport *

# TODO typo
cdef extern from "<OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>" namespace "OpenMS":

    cdef cppclass PeakFileOptions:

        PeakFileOptions() nogil except +
        PeakFileOptions(PeakFileOptions) nogil except +

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

        void setMz32Bit(bool mz_32_bit) nogil except +
        bool getMz32Bit() nogil except +
        void setIntensity32Bit(bool int_32_bit) nogil except +
        bool getIntensity32Bit() nogil except +

        void setRTRange(DRange1 & range_)
        bool hasRTRange()
        DRange1 getRTRange()
        # void setMZRange(DRange[ 1 ] & range_)
        # bool hasMZRange()
        # DRange[ 1 ]  getMZRange()
        # void setIntensityRange(DRange[ 1 ] & range_)
        # bool hasIntensityRange()
        # DRange[ 1 ]  getIntensityRange()

