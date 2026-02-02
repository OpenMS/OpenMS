from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from ExperimentalSettings cimport *
from smart_ptr cimport shared_ptr
from SwathMap cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/MRMFile.h>" namespace "OpenMS":

    cdef cppclass MRMFile:
        # wrap-doc:
        #  Minimal SRM/MRM file loader returning a single SwathMap wrapping the chromatogram container.
        MRMFile() except + nogil
        MRMFile(MRMFile &) except + nogil  # compiler

        libcpp_vector[ SwathMap ] loadMzML(String file_,
                                           String tmp,
                                           shared_ptr[ ExperimentalSettings ] exp_meta) except + nogil
