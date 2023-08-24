from Types cimport *
from libcpp cimport bool
from Types cimport *
from String cimport *
from InterfaceDataStructures cimport *
# from ISpectrumAccess cimport *
from IndexedMzMLDecoder cimport *
from MSSpectrum cimport *
from MSChromatogram cimport *

cdef extern from "<OpenMS/FORMAT/HANDLERS/IndexedMzMLHandler.h>" namespace "OpenMS":
    
    cdef cppclass IndexedMzMLHandler "OpenMS::Internal::IndexedMzMLHandler":
        IndexedMzMLHandler() except + nogil 
        IndexedMzMLHandler(IndexedMzMLHandler &) except + nogil 
        IndexedMzMLHandler(String filename) except + nogil 

        void openFile(String filename) except + nogil 
        bool getParsingSuccess() except + nogil 

        size_t getNrSpectra() except + nogil 
        size_t getNrChromatograms() except + nogil 

        shared_ptr[Spectrum] getSpectrumById(int id_) except + nogil 
        shared_ptr[Chromatogram] getChromatogramById(int id_) except + nogil 

        MSSpectrum getMSSpectrumById(int id_) except + nogil 
        void getMSSpectrumByNativeId(libcpp_string id_, MSSpectrum& spec) except + nogil 
        MSChromatogram getMSChromatogramById(int id_) except + nogil 
        void getMSChromatogramByNativeId(libcpp_string id_, MSChromatogram& chrom) except + nogil 

        void setSkipXMLChecks(bool skip) except + nogil 

