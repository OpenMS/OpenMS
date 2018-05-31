from String cimport *
from ExperimentalDesign cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/FORMAT/ExperimentalDesignFile.h>" namespace "OpenMS":
    cdef cppclass ExperimentalDesignFile "OpenMS::ExperimentalDesignFile":
        ExperimentalDesignFile() nogil except + #wrap-ignore
        ExperimentalDesignFile(ExperimentalDesignFile&) nogil except + #wrap-ignore
                
# COMMENT: wrap static methods
cdef extern from "<OpenMS/FORMAT/ExperimentalDesignFile.h>" namespace "OpenMS::ExperimentalDesignFile":
    ExperimentalDesign load(const String& tsv_file, bool) nogil except + #wrap-attach:ExperimentalDesignFile

