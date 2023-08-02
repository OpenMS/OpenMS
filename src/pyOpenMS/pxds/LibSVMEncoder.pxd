from Types cimport *
from ProgressLogger cimport *
from String cimport *
from SVMWrapper cimport *

cdef extern from "<OpenMS/FORMAT/LibSVMEncoder.h>" namespace "OpenMS":
    
    cdef cppclass LibSVMEncoder "OpenMS::LibSVMEncoder":
        # wrap-doc:
                #  Serves for encoding sequences into feature vectors
                #  
                #  The class can be used to construct composition vectors for
                #  sequences. Additionally the vectors can be encoded into
                #  the libsvm format

        LibSVMEncoder() except + nogil 
        LibSVMEncoder(LibSVMEncoder &) except + nogil  # compiler

# COMMENT: wrap static methods
cdef extern from "<OpenMS/FORMAT/LibSVMEncoder.h>" namespace "OpenMS::LibSVMEncoder":
        
        # static members
        libcpp_vector[double] predictPeptideRT(libcpp_vector[String] sequences,
                                               const SVMWrapper & svm,
                                               const String & allowed_characters,
                                               UInt maximum_sequence_length) except + nogil  # wrap-attach:LibSVMEncoder
