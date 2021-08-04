from Types cimport *
from libcpp cimport bool
from Map cimport *
from String cimport *

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecomposition.h>" namespace "OpenMS":
    
    cdef cppclass MassDecomposition "OpenMS::MassDecomposition":
        # wrap-doc:
            #   Class represents a decomposition of a mass into amino acids
            #   -----
            #   This class represents a mass decomposition into amino acids. A
            #   decomposition are amino acids given with frequencies which add
            #   up to a specific mass.

        MassDecomposition() nogil except +
        MassDecomposition(MassDecomposition &) nogil except +
        MassDecomposition(const String & deco) nogil except +
        # MassDecomposition  operator+=(MassDecomposition & d) nogil except +
        String toString() nogil except + # wrap-doc:Returns the decomposition as a string
        String toExpandedString() nogil except + # wrap-doc:Returns the decomposition as a string; instead of frequencies the amino acids are repeated
        # MassDecomposition operator+(MassDecomposition & rhs) nogil except +
        Size getNumberOfMaxAA() nogil except + # wrap-doc:Returns the max frequency of this composition
        # bool operator<(MassDecomposition & rhs) nogil except +
        # bool operator==(const String & deco) nogil except +
        bool containsTag(const String & tag) nogil except + # wrap-doc:Returns true if tag is contained in the mass decomposition
        bool compatible(MassDecomposition & deco) nogil except + # wrap-doc:Returns true if the mass decomposition if contained in this instance

