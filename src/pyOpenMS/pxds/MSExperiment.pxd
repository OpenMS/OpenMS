from libcpp.vector cimport vector as libcpp_vector
from MSSpectrum cimport *
from MSChromatogram cimport *
from DataValue cimport *
from String cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from ExperimentalSettings cimport *
from DateTime cimport *
from RangeManager cimport *

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/KERNEL/MSExperiment.h>" namespace "OpenMS":

    cdef cppclass MSExperiment(ExperimentalSettings, RangeManagerRtMzInt):
        # wrap-inherits:
        #   ExperimentalSettings
        #   RangeManagerRtMzInt
        #
        # wrap-doc:
        #   In-Memory representation of a mass spectrometry experiment.
        #   -----
        #   Contains the data and metadata of an experiment performed with an MS (or
        #   HPLC and MS). This representation of an MS experiment is organized as list
        #   of spectra and chromatograms and provides an in-memory representation of
        #   popular mass-spectrometric file formats such as mzXML or mzML. The
        #   meta-data associated with an experiment is contained in
        #   ExperimentalSettings (by inheritance) while the raw data (as well as
        #   spectra and chromatogram level meta data) is stored in objects of type
        #   MSSpectrum and MSChromatogram, which are accessible through the getSpectrum
        #   and getChromatogram functions.
        #   -----
        #   Spectra can be accessed by direct iteration or by getSpectrum(),
        #   while chromatograms are accessed through getChromatogram().
        #   See help(ExperimentalSettings) for information about meta-data.
        #   -----
        #   Usage:
        #     exp = MSExperiment()
        #     MzMLFile().load(path_to_file, exp)
        #     for spectrum in exp:
        #       print(spectrum.size()) # prints number of peaks
        #       mz, intensities = spectrum.get_peaks()
        #   -----

        MSExperiment() nogil except +
        MSExperiment(MSExperiment &) nogil except +

        ExperimentalSettings getExperimentalSettings() nogil except +
        
        # COMMENT: Spectra functions
        MSSpectrum& operator[](int) nogil except + # wrap-upper-limit:size()
        MSSpectrum getSpectrum(Size id_) nogil except + # wrap-ignore
        void addSpectrum(MSSpectrum spec) nogil except +
        void setSpectra(libcpp_vector[ MSSpectrum ] & spectra) nogil except +
        libcpp_vector[MSSpectrum] getSpectra() nogil except +
        void get2DPeakData(double min_rt, double max_rt, double min_mz, double max_mz, libcpp_vector[float] & rt, libcpp_vector[float] & mz, libcpp_vector[float] & intensity) nogil except + # wrap-ignore

        # COMMENT: Chromatogram functions
        MSChromatogram getChromatogram(Size id_) nogil except + # wrap-ignore
        void addChromatogram(MSChromatogram chromatogram) nogil except +
        void setChromatograms(libcpp_vector[MSChromatogram] chromatograms) nogil except +
        libcpp_vector[MSChromatogram] getChromatograms() nogil except +

        # COMMENT: Spectra iteration
        libcpp_vector[MSSpectrum].iterator begin() nogil except +        # wrap-iter-begin:__iter__(MSSpectrum)
        libcpp_vector[MSSpectrum].iterator end() nogil except +       # wrap-iter-end:__iter__(MSSpectrum)

        MSChromatogram calculateTIC() nogil except + # wrap-doc:Returns the total ion chromatogram
        void clear(bool clear_meta_data) nogil except + # wrap-doc:Clear all spectra data and meta data (if called with True)

        void updateRanges() nogil except + # wrap-doc:Recalculate global RT and m/z ranges after changes to the data has been made.
        void updateRanges(int msLevel) nogil except + # wrap-doc:Recalculate RT and m/z ranges for a specific MS level

        void reserveSpaceSpectra(Size s) nogil except +
        void reserveSpaceChromatograms(Size s) nogil except +

        # Size of experiment
        UInt64 getSize() nogil except + # wrap-doc:Returns the total number of peaks
        int size() nogil except +
        void resize(Size s) nogil except +
        bool empty() nogil except +
        void reserve(Size s) nogil except +
        Size getNrSpectra() nogil except + # wrap-doc:Returns the number of MS spectra
        Size getNrChromatograms() nogil except + # wrap-doc:Returns the number of chromatograms
        libcpp_vector[unsigned int] getMSLevels() nogil except +  # wrap-ignore

        void sortSpectra(bool sort_mz) nogil except + # wrap-doc:Sorts spectra by RT. If sort_mz=True also sort each peak in a spectrum by m/z
        void sortSpectra() nogil except + 
        void sortChromatograms(bool sort_rt) nogil except + # wrap-doc:Sorts chromatograms by m/z. If sort_rt=True also sort each chromatogram RT
        void sortChromatograms() nogil except +

        bool isSorted(bool check_mz) nogil except + # wrap-doc:Checks if all spectra are sorted with respect to ascending RT
        bool isSorted() nogil except +

        void getPrimaryMSRunPath(StringList& toFill) nogil except + # wrap-doc:References to the first MS file(s) after conversions. Used to trace results back to original data.
        void swap(MSExperiment) nogil except +

        bool operator==(MSExperiment) nogil except +
        void reset() nogil except +
        bool clearMetaDataArrays() nogil except +

        int getPrecursorSpectrum(int zero_based_index) nogil except + # wrap-doc:Returns the index of the precursor spectrum for spectrum at index @p zero_based_index
        
