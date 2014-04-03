from libcpp cimport bool

from CVMappings cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/CVMappingFile.h>" namespace "OpenMS":

    cdef cppclass CVMappingFile:

        CVMappingFile() nogil except +
        # CVMappingFile(CVMappingFile) nogil except + # private
        void load(String & filename, CVMappings & cv_mappings, bool strip_namespaces) nogil except +

