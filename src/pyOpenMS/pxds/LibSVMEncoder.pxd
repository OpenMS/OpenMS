from Types cimport *
from ProgressLogger cimport *
from String cimport *
from SVMWrapper cimport *

cdef extern from "<OpenMS/FORMAT/LibSVMEncoder.h>" namespace "OpenMS":
    
    cdef cppclass LibSVMEncoder "OpenMS::LibSVMEncoder":
        LibSVMEncoder() nogil except +
        LibSVMEncoder(LibSVMEncoder) nogil except + #wrap-ignore


# COMMENT: wrap static methods
cdef extern from "<OpenMS/FORMAT/LibSVMEncoder.h>" namespace "OpenMS::LibSVMEncoder":
        
        # static members
        libcpp_vector[double] predictPeptideRT(libcpp_vector[String] sequences,
                                               const String & allowed_characters,
                                               const SVMWrapper & svm,
                                               UInt maximum_sequence_length) nogil except + # wrap-attach:LibSVMEncoder
