from libcpp.vector cimport vector as libcpp_vector
from PeptideIdentification cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/DENOVO/CompNovoIdentificationCID.h>" namespace "OpenMS":

    cdef cppclass CompNovoIdentificationCID(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        CompNovoIdentificationCID()      nogil except +
        CompNovoIdentificationCID(CompNovoIdentificationCID)      nogil except + #private

        void getIdentifications(libcpp_vector[PeptideIdentification] & ids, MSExperiment[Peak1D,ChromatogramPeak])      nogil except +
        void getIdentification(PeptideIdentification & id, MSSpectrum[Peak1D] cid_spec)      nogil except +

