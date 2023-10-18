from Types cimport *
from libcpp cimport bool
from PeptideHit cimport *
from MSExperiment cimport *
from ResidueModification cimport *
from FASTAFile cimport *
from ProteaseDigestion cimport *

# cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS::OPXLDataStructs":
#    
#    cdef cppclass XLPrecursorComparator "OpenMS::OPXLDataStructs::XLPrecursorComparator":
#        XLPrecursorComparator(XLPrecursorComparator &) except + nogil 
#        bool operator()(XLPrecursor & a, XLPrecursor & b) except + nogil 
#        bool operator()(XLPrecursor & a, double & b) except + nogil 
#        bool operator()(double & a, XLPrecursor & b) except + nogil 
# 
