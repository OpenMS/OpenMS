from String cimport *
from ExperimentalDesign cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/FORMAT/ExperimentalDesignFile.h>" namespace "OpenMS":
    cdef cppclass ExperimentalDesignFile "OpenMS::ExperimentalDesignFile":
        ExperimentalDesignFile() except + nogil  # compiler
        ExperimentalDesignFile(ExperimentalDesignFile &) except + nogil  # compiler
                
# COMMENT: wrap static methods
cdef extern from "<OpenMS/FORMAT/ExperimentalDesignFile.h>" namespace "OpenMS::ExperimentalDesignFile":
    ExperimentalDesign load(const String& tsv_file, bool) except + nogil  #wrap-attach:ExperimentalDesignFile

