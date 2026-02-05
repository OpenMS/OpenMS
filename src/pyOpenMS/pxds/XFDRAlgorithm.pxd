from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from PeptideIdentification cimport *
from ProteinIdentification cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/XFDRAlgorithm.h>" namespace "OpenMS":

    cdef cppclass XFDRAlgorithm(DefaultParamHandler) :
        # wrap-doc:
        #  Calculates false discovery rate estimates on crosslink
        #  identifications
        # wrap-inherits:
        #  DefaultParamHandler
        XFDRAlgorithm() except + nogil
        XFDRAlgorithm(XFDRAlgorithm &) except + nogil

        ExitCodes run(PeptideIdentificationList& peptide_ids,
                      ProteinIdentification& protein_id) except + nogil

        ExitCodes validateClassArguments() except + nogil

cdef extern from "<OpenMS/ANALYSIS/XLMS/XFDRAlgorithm.h>" namespace "OpenMS::XFDRAlgorithm":
    cdef enum class ExitCodes "OpenMS::XFDRAlgorithm::ExitCodes":
        # wrap-attach:
        #    XFDRAlgorithm
        EXECUTION_OK
        ILLEGAL_PARAMETERS
        UNEXPECTED_RESULT
