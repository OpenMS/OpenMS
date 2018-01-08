from AbsoluteQuantitationStandards cimport *
from String cimport *
from Types cimport *

cdef extern from "<OpenMS/FORMAT/AbsoluteQuantitationStandardsFile.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitationStandardsFile:

        AbsoluteQuantitationStandardsFile() nogil except +
        AbsoluteQuantitationStandardsFile(AbsoluteQuantitationStandardsFile &) nogil except +

        void load(const String& filename, libcpp_vector[AQS_runConcentration]& run_concentrations) nogil except +
