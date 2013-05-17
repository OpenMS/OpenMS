from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from String cimport *

cdef extern from "<OpenMS/CHEMISTRY/PepIterator.h>" namespace "OpenMS":
    
    cdef cppclass PepIterator "OpenMS::PepIterator":
        # wrap-ignore
        # ABSTRACT class
        PepIterator() nogil except +
        PepIterator(PepIterator) nogil except +
        # POINTER # FASTAEntry operator*() nogil except +
        # PepIterator  operator++() nogil except +
        # POINTER # PepIterator * operator++(int ) nogil except +
        void setFastaFile(String & f) nogil except +
        String getFastaFile() nogil except +
        void setSpectrum(libcpp_vector[ double ] & s) nogil except +
        libcpp_vector[ double ]  getSpectrum() nogil except +
        void setTolerance(DoubleReal t) nogil except +
        DoubleReal getTolerance() nogil except +
        bool begin() nogil except +
        bool isAtEnd() nogil except +
        void registerChildren() nogil except +

