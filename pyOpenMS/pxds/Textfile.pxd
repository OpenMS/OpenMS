from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from StringList cimport *

cdef extern from "<OpenMS/FORMAT/TextFile.h>" namespace "OpenMS":
    
    cdef cppclass TextFile "OpenMS::TextFile":
        TextFile() nogil except +
        TextFile(TextFile) nogil except + #wrap-ignore
        TextFile(String &filename, bool trim_linesalse, Int first_n1) nogil except +
        void load(String &filename, bool trim_linesalse, Int first_n1) nogil except +
        void store(String &filename) nogil except +
