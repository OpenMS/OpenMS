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

        CompNovoIdentificationCID() nogil except +
        CompNovoIdentificationCID(CompNovoIdentificationCID &) nogil except +

        void getIdentifications(libcpp_vector[PeptideIdentification] & ids, MSExperiment) nogil except + # wrap-doc:Performs and returns de novo identifications
        void getIdentification(PeptideIdentification & id, MSSpectrum cid_spec) nogil except + # wrap-doc:Performs and returns de novo identifications
