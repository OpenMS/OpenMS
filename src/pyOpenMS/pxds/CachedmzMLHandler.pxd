from MSExperiment  cimport *
from MSSpectrum  cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from String cimport *
from ProgressLogger cimport *
from streampos cimport *

cdef extern from "<OpenMS/FORMAT/HANDLERS/CachedMzMLHandler.h>" namespace "OpenMS::Internal":

    # Do not use this class directly, rather use SpectrumAccessOpenMSCached

    cdef cppclass CachedMzMLHandler(ProgressLogger):
        # wrap-inherits:
        #   ProgressLogger

        CachedMzMLHandler() nogil except +
        CachedMzMLHandler(CachedMzMLHandler) nogil except +

        void writeMemdump(MSExperiment exp, String out) nogil except +
        void writeMetadata(MSExperiment exp, String out_meta) nogil except +

        void readMemdump(MSExperiment exp, String filename) nogil except +

        # void readSingleSpectrum(MSSpectrum & spectrum, String & filename, Size & idx) nogil except +
        libcpp_vector[ streampos ]  getSpectraIndex() nogil except +
        libcpp_vector[ streampos ]  getChromatogramIndex() nogil except +
        void createMemdumpIndex(String filename) nogil except +
        # NAMESPACE # void readSingleSpectrum(MSSpectrum & spectrum, std::ifstream & ifs, Size & idx)
        # NAMESPACE # void readSpectrumFast(OpenSwath::BinaryDataArrayPtr data1, OpenSwath::BinaryDataArrayPtr data2, std::ifstream & ifs, int ms_level, double rt)
        # NAMESPACE # void readChromatogramFast(OpenSwath::BinaryDataArrayPtr data1, OpenSwath::BinaryDataArrayPtr data2, std::ifstream & ifs)

