from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from ExperimentalSettings cimport *
from smart_ptr cimport shared_ptr
from SwathMap cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/TargetedDataFileLoader.h>" namespace "OpenMS":

    cdef cppclass TargetedDataFileLoader:
        # wrap-doc:
        #  Dispatcher that detects whether an mzML contains spectra (SWATH/DIA)
        #  or chromatograms only (SRM/MRM) and forwards to the appropriate loader.
        TargetedDataFileLoader() except + nogil
        TargetedDataFileLoader(TargetedDataFileLoader &) except + nogil  # compiler

        libcpp_vector[ SwathMap ] loadFile(String file_,
                                           String tmp,
                                           shared_ptr[ ExperimentalSettings ] exp_meta,
                                           String readoptions) except + nogil
