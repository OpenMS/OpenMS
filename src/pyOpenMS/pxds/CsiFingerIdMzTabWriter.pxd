from Types cimport *
from String cimport *
from MzTab cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>" namespace "OpenMS":
    
    cdef cppclass CsiFingerIdMzTabWriter "OpenMS::CsiFingerIdMzTabWriter":
        CsiFingerIdMzTabWriter() except + nogil  
        CsiFingerIdMzTabWriter(CsiFingerIdMzTabWriter &) except + nogil  # compiler

# wrap static method:
cdef extern from "<OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>" namespace "OpenMS::CsiFingerIdMzTabWriter":

        void read(libcpp_vector[ String ]& sirius_output_paths, 
                  const String& original_input_mzml, 
                  Size top_n_hits, 
                  MzTab& result) except + nogil  # wrap-attach:CsiFingerIdMzTabWriter

