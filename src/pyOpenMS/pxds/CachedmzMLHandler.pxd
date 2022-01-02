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

        CachedMzMLHandler() nogil except + # wrap-doc:An internal class that handles single spectra and chromatograms
        CachedMzMLHandler(CachedMzMLHandler &) nogil except + # compiler


        void writeMemdump(MSExperiment exp, String out) nogil except + # wrap-doc:Write complete spectra as a dump to the disk
        void writeMetadata(MSExperiment exp, String out_meta) nogil except + # wrap-doc:Write only the meta data of an MSExperiment

        void readMemdump(MSExperiment exp, String filename) nogil except + # wrap-doc:Read all spectra from a dump from the disk

        # void readSingleSpectrum(MSSpectrum & spectrum, String & filename, Size & idx) nogil except +
        libcpp_vector[ streampos ]  getSpectraIndex() nogil except +
        libcpp_vector[ streampos ]  getChromatogramIndex() nogil except +
        void createMemdumpIndex(String filename) nogil except + # wrap-doc:Create an index on the location of all the spectra and chromatograms
        # NAMESPACE # void readSingleSpectrum(MSSpectrum & spectrum, std::ifstream & ifs, Size & idx)
        # NAMESPACE # void readSpectrumFast(OpenSwath::BinaryDataArrayPtr data1, OpenSwath::BinaryDataArrayPtr data2, std::ifstream & ifs, int ms_level, double rt)
        # NAMESPACE # void readChromatogramFast(OpenSwath::BinaryDataArrayPtr data1, OpenSwath::BinaryDataArrayPtr data2, std::ifstream & ifs)
