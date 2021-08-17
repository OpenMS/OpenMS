from libcpp cimport bool

from CVMappings cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/CVMappingFile.h>" namespace "OpenMS":

    cdef cppclass CVMappingFile:

        CVMappingFile() nogil except +
        # private
        CVMappingFile(CVMappingFile) nogil except + # wrap-ignore

        void load(const String & filename, CVMappings & cv_mappings, bool strip_namespaces) nogil except + # wrap-doc:Loads CvMappings from the given file
