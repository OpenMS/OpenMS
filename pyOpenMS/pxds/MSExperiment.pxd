from libcpp.vector cimport vector as libcpp_vector
from MSSpectrum cimport *
from MSChromatogram cimport *
from DataValue cimport *
from String cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from MetaInfoInterface cimport *
from ExperimentalSettings cimport *
from DateTime cimport *

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/KERNEL/MSExperiment.h>" namespace "OpenMS":

    cdef cppclass MSExperiment[PeakT, ChromoPeakT](ExperimentalSettings):
        # wrap-inherits:
        #   ExperimentalSettings
        #
        # wrap-instances:
        #   MSExperiment := MSExperiment[Peak1D, ChromatogramPeak]

        MSExperiment() nogil except +
        MSExperiment(MSExperiment[PeakT, ChromoPeakT] &)  nogil except +

        bool operator==(MSExperiment[PeakT, ChromatogramPeak]) nogil except +
        void reset() nogil except +
        bool clearMetaDataArrays() nogil except +
        ExperimentalSettings getExperimentalSettings() nogil except +

        void swap(MSExperiment) nogil except +

        void setChromatograms(libcpp_vector[MSChromatogram[ChromoPeakT]] chromatograms) nogil except +
        void addChromatogram(MSChromatogram[ChromoPeakT] chromatogram) nogil except +
        libcpp_vector[MSChromatogram[ChromoPeakT]] getChromatograms() nogil except +

        MSChromatogram[ChromoPeakT] getTIC() nogil except +
        void clear(bool clear_meta_data)

        void   updateRanges() nogil except +
        void   updateRanges(int msLevel) nogil except +


        double getMinMZ() nogil except +
        double getMaxMZ() nogil except +
        double getMinRT() nogil except +
        double getMaxRT() nogil except +

        UInt64 getSize() nogil except +
        libcpp_vector[unsigned int] getMSLevels() nogil except +  # wrap-ignore


        void sortSpectra(bool sort_mz) nogil except +
        void sortSpectra() nogil except +
        void sortChromatograms(bool sort_rt) nogil except +
        void sortChromatograms() nogil except +

        bool isSorted(bool check_mz) nogil except +
        bool isSorted() nogil except +

        # from inheriting libcpp_vector[PeakT]:
        int   size() nogil except +
        MSSpectrum[PeakT] operator[](int)      nogil except + # wrap-upper-limit:size()
        void push_back(MSSpectrum[PeakT] spec) nogil except +

        libcpp_vector[MSSpectrum[PeakT]].iterator begin() nogil except +        # wrap-iter-begin:__iter__(MSSpectrum)
        libcpp_vector[MSSpectrum[PeakT]].iterator end()    nogil except +       # wrap-iter-end:__iter__(MSSpectrum)
        void  erase(libcpp_vector[MSSpectrum[PeakT]].iterator) nogil except +   # wrap-ignore

        # from MetaInfoInterface:


        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except +
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +


