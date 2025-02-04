from Types cimport *
from libcpp cimport bool
from String cimport *

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecomposition.h>" namespace "OpenMS":
    
    cdef cppclass MassDecomposition "OpenMS::MassDecomposition":
        # wrap-doc:
            #  Class represents a decomposition of a mass into amino acids
            #  
            #  This class represents a mass decomposition into amino acids. A
            #  decomposition are amino acids given with frequencies which add
            #  up to a specific mass.

        MassDecomposition() except + nogil 
        MassDecomposition(MassDecomposition &) except + nogil 
        MassDecomposition(const String & deco) except + nogil 
        # MassDecomposition  operator+=(MassDecomposition & d) except + nogil 
        String toString() except + nogil  # wrap-doc:Returns the decomposition as a string
        String toExpandedString() except + nogil  # wrap-doc:Returns the decomposition as a string; instead of frequencies the amino acids are repeated
        # MassDecomposition operator+(MassDecomposition & rhs) except + nogil 
        Size getNumberOfMaxAA() except + nogil  # wrap-doc:Returns the max frequency of this composition
        # bool operator<(MassDecomposition & rhs) except + nogil 
        # bool operator==(const String & deco) except + nogil 
        bool containsTag(const String & tag) except + nogil  # wrap-doc:Returns true if tag is contained in the mass decomposition
        bool compatible(MassDecomposition & deco) except + nogil  # wrap-doc:Returns true if the mass decomposition if contained in this instance

