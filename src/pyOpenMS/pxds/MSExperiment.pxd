from libcpp.vector cimport vector as libcpp_vector
from MSSpectrum cimport *
from MSChromatogram cimport *
from DataValue cimport *
from String cimport *
from Peak1D cimport *
from libcpp.pair cimport pair as libcpp_pair
from ChromatogramPeak cimport *
from ExperimentalSettings cimport *
from DateTime cimport *
from RangeManager cimport *
from Matrix cimport *

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/KERNEL/MSExperiment.h>" namespace "OpenMS":

    cdef cppclass MSExperiment(ExperimentalSettings, RangeManagerRtMzInt):
        # wrap-inherits:
        #  ExperimentalSettings
        #  RangeManagerRtMzInt
        #
        # wrap-doc:
        #  In-Memory representation of a mass spectrometry experiment.
        #  
        #  Contains the data and metadata of an experiment performed with an MS (or
        #  HPLC and MS). This representation of an MS experiment is organized as list
        #  of spectra and chromatograms and provides an in-memory representation of
        #  popular mass-spectrometric file formats such as mzXML or mzML. The
        #  meta-data associated with an experiment is contained in
        #  ExperimentalSettings (by inheritance) while the raw data (as well as
        #  spectra and chromatogram level meta data) is stored in objects of type
        #  MSSpectrum and MSChromatogram, which are accessible through the getSpectrum
        #  and getChromatogram functions.
        #  
        #  Spectra can be accessed by direct iteration or by getSpectrum(),
        #  while chromatograms are accessed through getChromatogram().
        #  See help(ExperimentalSettings) for information about meta-data.
        #  
        #  Usage:
        #
        #  .. code-block:: python
        #  
        #    exp = MSExperiment()
        #    MzMLFile().load(path_to_file, exp)
        #    for spectrum in exp:
        #      print(spectrum.size()) # prints number of peaks
        #      mz, intensities = spectrum.get_peaks()
        #  

        MSExperiment() except + nogil 
        MSExperiment(MSExperiment &) except + nogil 

        ExperimentalSettings getExperimentalSettings() except + nogil 
        
        # COMMENT: Spectra functions
        MSSpectrum& operator[](size_t) except + nogil  # wrap-upper-limit:size()
        MSSpectrum getSpectrum(Size id_) except + nogil  # wrap-ignore
        void addSpectrum(MSSpectrum spec) except + nogil 
        void setSpectra(libcpp_vector[ MSSpectrum ] & spectra) except + nogil 
        libcpp_vector[MSSpectrum] getSpectra() except + nogil 
        void get2DPeakData(double min_rt, double max_rt, double min_mz, double max_mz, unsigned int ms_level, libcpp_vector[float] & rt, libcpp_vector[float] & mz, libcpp_vector[float] & intensity) except + nogil  # wrap-ignore
        void get2DPeakDataIM(double min_rt, double max_rt, double min_mz, double max_mz, unsigned int ms_level, libcpp_vector[float] & rt, libcpp_vector[float] & mz, libcpp_vector[float] & intensity, libcpp_vector[float] & ion_mobility) except + nogil  # wrap-ignore
        void get2DPeakDataPerSpectrum(double min_rt, double max_rt, double min_mz, double max_mz, unsigned int ms_level, libcpp_vector[float] & rt, libcpp_vector[libcpp_vector[float]] & mz, libcpp_vector[libcpp_vector[float]] & intensity) except + nogil  # wrap-ignore
        void get2DPeakDataIMPerSpectrum(double min_rt, double max_rt, double min_mz, double max_mz, unsigned int ms_level, libcpp_vector[float] & rt, libcpp_vector[libcpp_vector[float]] & mz, libcpp_vector[libcpp_vector[float]] & intensity, libcpp_vector[libcpp_vector[float]] & ion_mobility) except + nogil  # wrap-ignore
        libcpp_vector[libcpp_vector[double]] aggregateFromMatrix(Matrix[double] & ranges, unsigned int ms_level, libcpp_string mz_agg) except + nogil # wrap-doc:Aggregates intensity values for multiple m/z and RT ranges specified in a matrix
        libcpp_vector[MSChromatogram] extractXICsFromMatrix(Matrix[double] & ranges, unsigned int ms_level, libcpp_string mz_agg) except + nogil # wrap-doc:Extracts XIC chromatograms for multiple m/z and RT ranges specified in a matrix

        # COMMENT: Chromatogram functions
        MSChromatogram getChromatogram(Size id_) except + nogil  # wrap-ignore
        void addChromatogram(MSChromatogram chromatogram) except + nogil 
        void setChromatograms(libcpp_vector[MSChromatogram] chromatograms) except + nogil 
        libcpp_vector[MSChromatogram] getChromatograms() except + nogil 

        # COMMENT: Spectra iteration
        libcpp_vector[MSSpectrum].iterator begin() except + nogil         # wrap-iter-begin:__iter__(MSSpectrum)
        libcpp_vector[MSSpectrum].iterator end() except + nogil        # wrap-iter-end:__iter__(MSSpectrum)

        MSChromatogram calculateTIC() except + nogil  # wrap-doc:Returns the total ion chromatogram
        void clear(bool clear_meta_data) except + nogil  # wrap-doc:Clear all spectra data and meta data (if called with True)

        void updateRanges() except + nogil  # wrap-doc:Recalculate global RT and m/z ranges after changes to the data has been made.
        void updateRanges(int msLevel) except + nogil  # wrap-doc:Recalculate RT and m/z ranges for a specific MS level

        void reserveSpaceSpectra(Size s) except + nogil 
        void reserveSpaceChromatograms(Size s) except + nogil 

        # Size of experiment
        UInt64 getSize() except + nogil  # wrap-doc:Returns the total number of peaks
        int size() except + nogil 
        void resize(Size s) except + nogil 
        bool empty() except + nogil 
        void reserve(Size s) except + nogil 
        Size getNrSpectra() except + nogil  # wrap-doc:Returns the number of MS spectra
        Size getNrChromatograms() except + nogil  # wrap-doc:Returns the number of chromatograms
        libcpp_vector[unsigned int] getMSLevels() except + nogil   # wrap-ignore

        void sortSpectra(bool sort_mz) except + nogil  # wrap-doc:Sorts spectra by RT. If sort_mz=True also sort each peak in a spectrum by m/z
        void sortSpectra() except + nogil  
        void sortChromatograms(bool sort_rt) except + nogil  # wrap-doc:Sorts chromatograms by m/z. If sort_rt=True also sort each chromatogram RT
        void sortChromatograms() except + nogil 

        bool isSorted(bool check_mz) except + nogil  # wrap-doc:Checks if all spectra are sorted with respect to ascending RT
        bool isSorted() except + nogil 

        void getPrimaryMSRunPath(StringList& toFill) except + nogil  # wrap-doc:References to the first MS file(s) after conversions. Used to trace results back to original data.
        void swap(MSExperiment) except + nogil 

        bool operator==(MSExperiment) except + nogil 
        void reset() except + nogil 
        bool clearMetaDataArrays() except + nogil 

        int getPrecursorSpectrum(int zero_based_index) except + nogil  # wrap-doc:Returns the index of the precursor spectrum for spectrum at index @p zero_based_index
        
