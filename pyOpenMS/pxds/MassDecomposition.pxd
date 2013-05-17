from Types cimport *
from libcpp cimport bool
from Map cimport *
from String cimport *

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecomposition.h>" namespace "OpenMS":
    
    cdef cppclass MassDecomposition "OpenMS::MassDecomposition":
        MassDecomposition() nogil except +
        MassDecomposition(MassDecomposition) nogil except +
        MassDecomposition(String & deco) nogil except +
        # MassDecomposition  operator+=(MassDecomposition & d) nogil except +
        String toString() nogil except +
        String toExpandedString() nogil except +
        # MassDecomposition operator+(MassDecomposition & rhs) nogil except +
        Size getNumberOfMaxAA() nogil except +
        # bool operator<(MassDecomposition & rhs) nogil except +
        # bool operator==(String & deco) nogil except +
        bool containsTag(String & tag) nogil except +
        bool compatible(MassDecomposition & deco) nogil except +

