from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from PeptideIdentification cimport *
from ProteinIdentification cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/XFDRAlgorithm.h>" namespace "OpenMS":

    cdef cppclass XFDRAlgorithm(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        XFDRAlgorithm() except + nogil 
        XFDRAlgorithm(XFDRAlgorithm &) except + nogil 

        XFDRAlgorithm_ExitCodes run(libcpp_vector[ PeptideIdentification ]& peptide_ids,
                                    ProteinIdentification& protein_id) except + nogil 

        XFDRAlgorithm_ExitCodes validateClassArguments() except + nogil 

cdef extern from "<OpenMS/ANALYSIS/XLMS/XFDRAlgorithm.h>" namespace "OpenMS::XFDRAlgorithm":
    cdef enum XFDRAlgorithm_ExitCodes "OpenMS::XFDRAlgorithm::ExitCodes":
        #wrap-attach:
        #   XFDRAlgorithm
        EXECUTION_OK
        ILLEGAL_PARAMETERS
        UNEXPECTED_RESULT
