from libcpp cimport bool

from CVMappings cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/CVMappingFile.h>" namespace "OpenMS":

    cdef cppclass CVMappingFile:

        CVMappingFile() except + nogil 
        # private
        CVMappingFile(CVMappingFile) except + nogil  # wrap-ignore

        void load(const String & filename, CVMappings & cv_mappings, bool strip_namespaces) except + nogil  # wrap-doc:Loads CvMappings from the given file
