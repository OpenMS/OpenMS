from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from DRange cimport *
from MSNumpressCoder cimport *

cdef extern from "<OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>" namespace "OpenMS":

    cdef cppclass PeakFileOptions:
        # wrap-doc:
        #  Options for loading files containing peak data
        PeakFileOptions() except + nogil 
        PeakFileOptions(PeakFileOptions &) except + nogil 

        void setMetadataOnly(bool) except + nogil  # wrap-doc:Sets whether or not to load only meta data
        bool getMetadataOnly()     except + nogil  # wrap-doc:Returns whether or not to load only meta data

        void setWriteSupplementalData(bool) except + nogil  # wrap-doc:Sets whether or not to write supplemental peak data in MzData files
        bool getWriteSupplementalData()     except + nogil  # wrap-doc:Returns whether or not to write supplemental peak data in MzData files

        void setMSLevels(libcpp_vector[int] levels) except + nogil  # wrap-doc:Sets the desired MS levels for peaks to load
        void addMSLevel(Int level) except + nogil  # wrap-doc:Adds a desired MS level for peaks to load
        void clearMSLevels()       except + nogil  # wrap-doc:Clears the MS levels
        bool hasMSLevels()         except + nogil  # wrap-doc:Returns true, if MS levels have been set
        bool containsMSLevel(int level)  except + nogil  # wrap-doc:Returns true, if MS level `level` has been set
        libcpp_vector[int] getMSLevels()    except + nogil  # wrap-doc:Returns the set MS levels

        void setCompression(bool) except + nogil  # wrap-doc:Sets if data should be compressed when writing
        bool getCompression()     except + nogil  # wrap-doc:Returns true, if data should be compressed when writing

        void setMz32Bit(bool mz_32_bit) except + nogil  # wrap-doc:Sets if mz-data and rt-data should be stored with 32bit or 64bit precision
        bool getMz32Bit() except + nogil  # wrap-doc:Returns true, if mz-data and rt-data should be stored with 32bit precision
        void setIntensity32Bit(bool int_32_bit) except + nogil  # wrap-doc:Sets if intensity data should be stored with 32bit or 64bit precision
        bool getIntensity32Bit() except + nogil  # wrap-doc:Returns true, if intensity data should be stored with 32bit precision

        void setRTRange(DRange1 & range_) except + nogil  # wrap-doc:Restricts the range of RT values for peaks to load
        bool hasRTRange() except + nogil  # wrap-doc:Returns true if an RT range has been set
        DRange1 getRTRange() except + nogil  # wrap-doc:Returns the RT range
        void setMZRange(DRange1 & range_) except + nogil  # wrap-doc:Restricts the range of MZ values for peaks to load
        bool hasMZRange() except + nogil  # wrap-doc:Returns true if an MZ range has been set
        DRange1 getMZRange() except + nogil  # wrap-doc:Returns the MZ range
        void setIntensityRange(DRange1 & range_) except + nogil  # wrap-doc:Restricts the range of intensity values for peaks to load
        bool hasIntensityRange() except + nogil  # wrap-doc:Returns true if an intensity range has been set
        DRange1 getIntensityRange() except + nogil  # wrap-doc:Returns the intensity range

        Size getMaxDataPoolSize() except + nogil  # wrap-doc:Returns maximal size of the data pool
        void setMaxDataPoolSize(Size s) except + nogil  # wrap-doc:Sets maximal size of the data pool

        void setSortSpectraByMZ(bool doSort) except + nogil  # wrap-doc:Sets whether or not to sort peaks in spectra
        bool getSortSpectraByMZ() except + nogil  # wrap-doc:Returns whether or not peaks in spectra should be sorted
        void setSortChromatogramsByRT(bool doSort) except + nogil  # wrap-doc:Sets whether or not to sort peaks in chromatograms
        bool getSortChromatogramsByRT() except + nogil  # wrap-doc:Returns whether or not peaks in chromatograms should be sorted

        bool hasFilters() except + nogil 

        void setFillData(bool only) except + nogil  # wrap-doc:Sets whether to fill the actual data into the container (spectrum/chromatogram)
        bool getFillData() except + nogil  # wrap-doc:Returns whether to fill the actual data into the container (spectrum/chromatogram)

        void setSkipXMLChecks(bool only) except + nogil  # wrap-doc:Sets whether to skip some XML checks and be fast instead
        bool getSkipXMLChecks() except + nogil  # wrap-doc:Returns whether to skip some XML checks and be fast instead

        bool getWriteIndex() except + nogil  # wrap-doc:Returns whether to write an index at the end of the file (e.g. indexedmzML file format)
        void setWriteIndex(bool write_index) except + nogil  # wrap-doc:Returns whether to write an index at the end of the file (e.g. indexedmzML file format)

        NumpressConfig getNumpressConfigurationMassTime() except + nogil  # wrap-doc:Sets numpress configuration options for m/z or rt dimension
        void setNumpressConfigurationMassTime(NumpressConfig config) except + nogil  # wrap-doc:Returns numpress configuration options for m/z or rt dimension

        NumpressConfig getNumpressConfigurationIntensity() except + nogil  # wrap-doc:Sets numpress configuration options for intensity dimension
        void setNumpressConfigurationIntensity(NumpressConfig config) except + nogil  # wrap-doc:Returns numpress configuration options for intensity dimension

        NumpressConfig getNumpressConfigurationFloatDataArray() except + nogil  # wrap-doc:Sets numpress configuration options for float data arrays
        void setNumpressConfigurationFloatDataArray(NumpressConfig config) except + nogil  # wrap-doc:Returns numpress configuration options for float data arrays

        void setForceMQCompatability(bool forceMQ) except + nogil  # wrap-doc:[mzXML only!]Returns Whether to write a scan-index and meta data to indicate a Thermo FTMS/ITMS instrument (required to have parameter control in MQ)
        bool getForceMQCompatability() except + nogil  # wrap-doc:[mzXML only!]Returns Whether to write a scan-index and meta data to indicate a Thermo FTMS/ITMS instrument (required to have parameter control in MQ)

        void setForceTPPCompatability(bool forceTPP) except + nogil  # wrap-doc:[ mzML only!]Returns Whether to skip writing the \<isolationWindow\> tag so that TPP finds the correct precursor m/z
        bool getForceTPPCompatability() except + nogil  # wrap-doc:[mzML only!]Returns Whether to skip writing the \<isolationWindow\> tag so that TPP finds the correct precursor m/z
        
