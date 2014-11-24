from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from PepIterator cimport *

cdef extern from "<OpenMS/CHEMISTRY/TrypticIterator.h>" namespace "OpenMS":
    
    cdef cppclass TrypticIterator(PepIterator) :
        # wrap-inherits:
        #  PepIterator

        TrypticIterator() nogil except +
        TrypticIterator(TrypticIterator) nogil except +
        # FASTAEntry operator*() nogil except +
        # PepIterator  operator++() nogil except +
        # POINTER # PepIterator * operator++(int i) nogil except +

        # in base class

        # void setFastaFile(String & f) nogil except +
        # String getFastaFile() nogil except +
        # void setTolerance(double ) nogil except +
        # double getTolerance() nogil except +
        # void setSpectrum(libcpp_vector[ double ] & ) nogil except +
        # libcpp_vector[ double ]  getSpectrum() nogil except +
        # bool begin() nogil except +
        # bool isAtEnd() nogil except +
        bool isDigestingEnd(char aa1, char aa2) nogil except +
        String getProductName() nogil except +
        # POINTER # PepIterator * create() nogil except +

