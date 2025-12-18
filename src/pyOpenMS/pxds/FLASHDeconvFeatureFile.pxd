from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from libcpp cimport bool
from Types cimport *
from String cimport *
from DeconvolvedSpectrum cimport *
from FLASHHelperClasses cimport *

cdef extern from "<OpenMS/FORMAT/FLASHDeconvFeatureFile.h>" namespace "OpenMS":

    cdef cppclass FLASHDeconvFeatureFile:
        # wrap-doc:
        #  FLASHDeconv feature level output *.tsv, *.ms1ft (for Promex), *.feature (for TopPIC) file formats.
        #  This class provides static methods for writing mass feature data.
        #  Note: Methods taking std::ostream are not directly exposed. Use file-based workflows.

        FLASHDeconvFeatureFile() except + nogil
        FLASHDeconvFeatureFile(FLASHDeconvFeatureFile &) except + nogil

