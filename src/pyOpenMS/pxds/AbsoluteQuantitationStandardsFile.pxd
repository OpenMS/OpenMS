from AbsoluteQuantitationStandards cimport *
from String cimport *
from Types cimport *

cdef extern from "<OpenMS/FORMAT/AbsoluteQuantitationStandardsFile.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitationStandardsFile:
        # wrap-doc:
        #  Load files containing runConcentration data

        AbsoluteQuantitationStandardsFile() except + nogil 
        AbsoluteQuantitationStandardsFile(AbsoluteQuantitationStandardsFile &) except + nogil  # compiler

        void load(const String& filename, libcpp_vector[AQS_runConcentration]& run_concentrations) except + nogil 
