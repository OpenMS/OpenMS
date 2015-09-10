from Types cimport *
from libcpp cimport bool
from Types cimport *
from String cimport *
from InterfaceDataStructures cimport *
# from ISpectrumAccess cimport *
from IndexedMzMLDecoder cimport *

cdef extern from "<OpenMS/FORMAT/IndexedMzMLFile.h>" namespace "OpenMS":
    
    cdef cppclass IndexedMzMLFile "OpenMS::IndexedMzMLFile":
        IndexedMzMLFile() nogil except +
        IndexedMzMLFile(IndexedMzMLFile) nogil except +
        IndexedMzMLFile(String filename) nogil except +
        void openFile(String filename) nogil except +
        bool getParsingSuccess() nogil except +
        size_t getNrSpectra() nogil except +
        size_t getNrChromatograms() nogil except +
        shared_ptr[Spectrum] getSpectrumById(int id_) nogil except +
        shared_ptr[Chromatogram] getChromatogramById(int id_) nogil except +

