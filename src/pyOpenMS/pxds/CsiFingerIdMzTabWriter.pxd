from Types cimport *
from String cimport *
from MzTab cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>" namespace "OpenMS":
    
    cdef cppclass CsiFingerIdMzTabWriter "OpenMS::CsiFingerIdMzTabWriter":
        CsiFingerIdMzTabWriter() nogil except + 
        CsiFingerIdMzTabWriter(CsiFingerIdMzTabWriter) nogil except +

# wrap static method:
cdef extern from "<OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>" namespace "OpenMS::CsiFingerIdMzTabWriter":

        void read(libcpp_vector[ String ]& sirius_output_paths, 
                  const String& original_input_mzml, 
                  Size top_n_hits, 
                  MzTab& result) nogil except + # wrap-attach:CsiFingerIdMzTabWriter

