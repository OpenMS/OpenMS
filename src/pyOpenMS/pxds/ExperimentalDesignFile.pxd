from String cimport *
from ExperimentalDesign cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/FORMAT/ExperimentalDesignFile.h>" namespace "OpenMS":
    cdef cppclass ExperimentalDesignFile "OpenMS::ExperimentalDesignFile":
        # wrap-doc:
        #  Load an experimental design from a TSV file. (see ExperimentalDesign
        #  for details on the supported format)
        ExperimentalDesignFile() except + nogil  # compiler
        ExperimentalDesignFile(ExperimentalDesignFile &) except + nogil  # compiler

        @staticmethod
        ExperimentalDesign load(const String& tsv_file, bool) except + nogil
