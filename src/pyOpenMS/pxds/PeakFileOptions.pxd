from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from DRange cimport *
from MSNumpressCoder cimport *

cdef extern from "<OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>" namespace "OpenMS":

    cdef cppclass PeakFileOptions:
        # wrap-doc:
        #   Options for loading files containing peak data
        PeakFileOptions() nogil except +
        PeakFileOptions(PeakFileOptions &) nogil except +

        void setMetadataOnly(bool) nogil except + # wrap-doc:Sets whether or not to load only meta data
        bool getMetadataOnly()     nogil except + # wrap-doc:Returns whether or not to load only meta data

        void setWriteSupplementalData(bool) nogil except + # wrap-doc:Sets whether or not to write supplemental peak data in MzData files
        bool getWriteSupplementalData()     nogil except + # wrap-doc:Returns whether or not to write supplemental peak data in MzData files

        void setMSLevels(libcpp_vector[int] levels) nogil except + # wrap-doc:Sets the desired MS levels for peaks to load
        void addMSLevel(Int level) nogil except + # wrap-doc:Adds a desired MS level for peaks to load
        void clearMSLevels()       nogil except + # wrap-doc:Clears the MS levels
        bool hasMSLevels()         nogil except + # wrap-doc:Returns true, if MS levels have been set
        bool containsMSLevel(int level)  nogil except + # wrap-doc:Returns true, if MS level `level` has been set
        libcpp_vector[int] getMSLevels()    nogil except + # wrap-doc:Returns the set MS levels

        void setCompression(bool) nogil except + # wrap-doc:Sets if data should be compressed when writing
        bool getCompression()     nogil except + # wrap-doc:Returns true, if data should be compressed when writing

        void setMz32Bit(bool mz_32_bit) nogil except + # wrap-doc:Sets if mz-data and rt-data should be stored with 32bit or 64bit precision
        bool getMz32Bit() nogil except + # wrap-doc:Returns true, if mz-data and rt-data should be stored with 32bit precision
        void setIntensity32Bit(bool int_32_bit) nogil except + # wrap-doc:Sets if intensity data should be stored with 32bit or 64bit precision
        bool getIntensity32Bit() nogil except + # wrap-doc:Returns true, if intensity data should be stored with 32bit precision

        void setRTRange(DRange1 & range_) nogil except + # wrap-doc:Restricts the range of RT values for peaks to load
        bool hasRTRange() nogil except + # wrap-doc:Returns true if an RT range has been set
        DRange1 getRTRange() nogil except + # wrap-doc:Returns the RT range
        void setMZRange(DRange1 & range_) nogil except + # wrap-doc:Restricts the range of MZ values for peaks to load
        bool hasMZRange() nogil except + # wrap-doc:Returns true if an MZ range has been set
        DRange1 getMZRange() nogil except + # wrap-doc:Returns the MZ range
        void setIntensityRange(DRange1 & range_) nogil except + # wrap-doc:Restricts the range of intensity values for peaks to load
        bool hasIntensityRange() nogil except + # wrap-doc:Returns true if an intensity range has been set
        DRange1 getIntensityRange() nogil except + # wrap-doc:Returns the intensity range

        Size getMaxDataPoolSize() nogil except + # wrap-doc:Returns maximal size of the data pool
        void setMaxDataPoolSize(Size s) nogil except + # wrap-doc:Sets maximal size of the data pool

        void setSortSpectraByMZ(bool doSort) nogil except + # wrap-doc:Sets whether or not to sort peaks in spectra
        bool getSortSpectraByMZ() nogil except + # wrap-doc:Returns whether or not peaks in spectra should be sorted
        void setSortChromatogramsByRT(bool doSort) nogil except + # wrap-doc:Sets whether or not to sort peaks in chromatograms
        bool getSortChromatogramsByRT() nogil except + # wrap-doc:Returns whether or not peaks in chromatograms should be sorted

        bool hasFilters() nogil except +

        void setFillData(bool only) nogil except + # wrap-doc:Sets whether to fill the actual data into the container (spectrum/chromatogram)
        bool getFillData() nogil except + # wrap-doc:Returns whether to fill the actual data into the container (spectrum/chromatogram)

        void setSkipXMLChecks(bool only) nogil except + # wrap-doc:Sets whether to skip some XML checks and be fast instead
        bool getSkipXMLChecks() nogil except + # wrap-doc:Returns whether to skip some XML checks and be fast instead

        bool getWriteIndex() nogil except + # wrap-doc:Returns whether to write an index at the end of the file (e.g. indexedmzML file format)
        void setWriteIndex(bool write_index) nogil except + # wrap-doc:Returns whether to write an index at the end of the file (e.g. indexedmzML file format)

        NumpressConfig getNumpressConfigurationMassTime() nogil except + # wrap-doc:Sets numpress configuration options for m/z or rt dimension
        void setNumpressConfigurationMassTime(NumpressConfig config) nogil except + # wrap-doc:Returns numpress configuration options for m/z or rt dimension

        NumpressConfig getNumpressConfigurationIntensity() nogil except + # wrap-doc:Sets numpress configuration options for intensity dimension
        void setNumpressConfigurationIntensity(NumpressConfig config) nogil except + # wrap-doc:Returns numpress configuration options for intensity dimension

        NumpressConfig getNumpressConfigurationFloatDataArray() nogil except + # wrap-doc:Sets numpress configuration options for float data arrays
        void setNumpressConfigurationFloatDataArray(NumpressConfig config) nogil except + # wrap-doc:Returns numpress configuration options for float data arrays

        void setForceMQCompatability(bool forceMQ) nogil except + # wrap-doc:[mzXML only!]Returns Whether to write a scan-index and meta data to indicate a Thermo FTMS/ITMS instrument (required to have parameter control in MQ)
        bool getForceMQCompatability() nogil except + # wrap-doc:[mzXML only!]Returns Whether to write a scan-index and meta data to indicate a Thermo FTMS/ITMS instrument (required to have parameter control in MQ)

        void setForceTPPCompatability(bool forceTPP) nogil except + # wrap-doc:[ mzML only!]Returns Whether to skip writing the \<isolationWindow\> tag so that TPP finds the correct precursor m/z
        bool getForceTPPCompatability() nogil except + # wrap-doc:[mzML only!]Returns Whether to skip writing the \<isolationWindow\> tag so that TPP finds the correct precursor m/z
        
