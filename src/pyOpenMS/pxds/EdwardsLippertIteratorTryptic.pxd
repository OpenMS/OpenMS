from Types cimport *
from libcpp cimport bool
from EdwardsLippertIterator cimport *

cdef extern from "<OpenMS/CHEMISTRY/EdwardsLippertIteratorTryptic.h>" namespace "OpenMS":
    
    cdef cppclass EdwardsLippertIteratorTryptic(EdwardsLippertIterator) :
        # wrap-inherits:
        #  EdwardsLippertIterator
        EdwardsLippertIteratorTryptic() nogil except +
        EdwardsLippertIteratorTryptic(EdwardsLippertIteratorTryptic) nogil except +

