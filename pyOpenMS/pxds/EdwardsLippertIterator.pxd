from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from PepIterator cimport *

cdef extern from "<OpenMS/CHEMISTRY/EdwardsLippertIterator.h>" namespace "OpenMS":
    
    cdef cppclass EdwardsLippertIterator(PepIterator) :
        # wrap-inherits:
        #  PepIterator
        EdwardsLippertIterator() nogil except +
        EdwardsLippertIterator(EdwardsLippertIterator) nogil except +
        bool isDigestingEnd(char AA1, char AA2) nogil except +
        String getProductName() nogil except +
        # POINTER # PepIterator * create() nogil except +

