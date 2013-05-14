from MSSpectrum cimport *
from MSChromatogram cimport *
from DataValue cimport *
from String cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from MetaInfoInterface cimport *

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/KERNEL/MSExperiment.h>" namespace "OpenMS":

    cdef cppclass MSExperiment[PeakT, ChromoPeakT](MetaInfoInterface):
        # wrap-inherits:
        #   MetaInfoInterface
        #
        # wrap-instances:
        #   MSExperiment := MSExperiment[Peak1D, ChromatogramPeak]

        MSExperiment() nogil except +
        MSExperiment(MSExperiment[PeakT, ChromoPeakT] &)  nogil except +

        double getMinMZ() nogil except +
        double getMaxMZ() nogil except +
        double getMinRT() nogil except +
        double getMaxRT() nogil except +

        void sortSpectra(bool) nogil except +
        bool isSorted(bool check_mz) nogil except +
        bool isSorted() nogil except +

        int   size() nogil except +

        MSSpectrum[PeakT] operator[](int)      nogil except + # wrap-upper-limit:size()

        void   updateRanges() nogil except +
        void   updateRanges(int msLevel) nogil except +

        UInt64 getSize() nogil except +


        void push_back(MSSpectrum[PeakT] spec) nogil except +
        String getLoadedFilePath() nogil except +
        void setLoadedFilePath(String path) nogil except +

        libcpp_vector[MSSpectrum[PeakT]].iterator begin() nogil except +        # wrap-iter-begin:__iter__(MSSpectrum)
        libcpp_vector[MSSpectrum[PeakT]].iterator end()    nogil except +       # wrap-iter-end:__iter__(MSSpectrum)
        void  erase(libcpp_vector[MSSpectrum[PeakT]].iterator) nogil except +   # wrap-ignore

        libcpp_vector[MSChromatogram[ChromatogramPeak]] getChromatograms() nogil except +
        void setChromatograms(libcpp_vector[MSChromatogram[ChromatogramPeak]] chromatograms) nogil except +
        void addChromatogram(MSChromatogram[ChromatogramPeak] chromatogram) nogil except +

        void clear(bool clear_meta_data)

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


