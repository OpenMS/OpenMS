from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from DRange cimport *
from MSNumpressCoder cimport *

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

        void setRTRange(DRange1 & range_) nogil except +
        bool hasRTRange() nogil except +
        DRange1 getRTRange() nogil except +
        void setMZRange(DRange1 & range_) nogil except +
        bool hasMZRange() nogil except +
        DRange1 getMZRange() nogil except +
        void setIntensityRange(DRange1 & range_) nogil except +
        bool hasIntensityRange() nogil except +
        DRange1 getIntensityRange() nogil except +

        Size getMaxDataPoolSize() nogil except +
        void setMaxDataPoolSize(Size s) nogil except +

        void setSortSpectraByMZ(bool doSort) nogil except +
        bool getSortSpectraByMZ() nogil except +
        void setSortChromatogramsByRT(bool doSort) nogil except +
        bool getSortChromatogramsByRT() nogil except +

        bool hasFilters() nogil except +

        void setFillData(bool only) nogil except +
        bool getFillData() nogil except +

        void setSkipXMLChecks(bool only) nogil except +
        bool getSkipXMLChecks() nogil except +

        bool getWriteIndex() nogil except +
        void setWriteIndex(bool write_index) nogil except +

        NumpressConfig getNumpressConfigurationMassTime() nogil except +
        void setNumpressConfigurationMassTime(NumpressConfig config) nogil except +

        NumpressConfig getNumpressConfigurationIntensity() nogil except +
        void setNumpressConfigurationIntensity(NumpressConfig config) nogil except +

        NumpressConfig getNumpressConfigurationFloatDataArray() nogil except +
        void setNumpressConfigurationFloatDataArray(NumpressConfig config) nogil except +

        void setForceMQCompatability(bool forceMQ) nogil except +
        bool getForceMQCompatability() nogil except +

        void setForceTPPCompatability(bool forceTPP) nogil except +
        bool getForceTPPCompatability() nogil except +
        
