from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from ProgressLogger cimport *
from String cimport *
from TextFile cimport *
from File cimport *

cdef extern from "<OpenMS/ANALYSIS/SVM/SVMWrapper.h>" namespace "OpenMS":
    
    cdef cppclass SVMData "OpenMS::SVMData":
        SVMData() nogil except +
        SVMData(SVMData) nogil except + #wrap-ignore
        # libcpp_vector[ libcpp_vector[ std::pair[ Int, DoubleReal ] ] ] sequences
        libcpp_vector[ double ] labels
        # SVMData(libcpp_vector[ libcpp_vector[ std::pair[ Int, DoubleReal ] ] ] &seqs, libcpp_vector[ DoubleReal ] &lbls) nogil except +
        bool operator==(SVMData &rhs) nogil except +
        bool store(String &filename) nogil except +
        bool load(String &filename) nogil except +

