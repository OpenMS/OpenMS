from Types cimport *
from AbsoluteQuantitationMethod cimport *

cdef extern from "<OpenMS/FORMAT/AbsoluteQuantitationMethodFile.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitationMethodFile:
        # wrap-doc:
        #  File adapter for AbsoluteQuantitationMethod files. Loads and stores .csv or .tsv
        #  files describing absolute quantitation methods including calibration curve
        #  parameters, limits of detection/quantitation, and transformation models

        AbsoluteQuantitationMethodFile() except + nogil
        AbsoluteQuantitationMethodFile(AbsoluteQuantitationMethodFile &) except + nogil  # compiler

        void load(const String& filename, libcpp_vector[ AbsoluteQuantitationMethod ]& aqm_list) except + nogil  # wrap-doc:Loads AbsoluteQuantitationMethod data from a CSV/TSV file into a list of methods
        void store(const String& filename, const libcpp_vector[ AbsoluteQuantitationMethod ]& aqm_list) except + nogil  # wrap-doc:Stores a list of AbsoluteQuantitationMethod objects to a CSV/TSV file
