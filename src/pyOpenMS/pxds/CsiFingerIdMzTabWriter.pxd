from Types cimport *
from String cimport *
from MzTab cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>" namespace "OpenMS":
    
    cdef cppclass CsiFingerIdMzTabWriter "OpenMS::CsiFingerIdMzTabWriter":
        CsiFingerIdMzTabWriter() nogil except +  # wrap-ignore
        CsiFingerIdMzTabWriter(CsiFingerIdMzTabWriter) nogil except + #wrap-ignore

#
# wrap static method:
#
cdef extern from "<OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>" namespace "OpenMS::CsiFingerIdMzTabWriter":

   void read(libcpp_vector[String] sirius_output_paths, String original_input_mzml, Size top_n_hits, MzTab & result) nogil except + # wrap-attach:CsiFingerIdMzTabWriter

