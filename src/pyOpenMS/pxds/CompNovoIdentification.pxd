from libcpp.vector cimport vector as libcpp_vector
from PeptideIdentification cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/DENOVO/CompNovoIdentification.h>" namespace "OpenMS":

    cdef cppclass CompNovoIdentification(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        CompNovoIdentification() nogil except +
        CompNovoIdentification(CompNovoIdentification &) nogil except +

        void getIdentifications(libcpp_vector[PeptideIdentification] & ids, MSExperiment) nogil except + # wrap-doc:Performs and returns de novo identifications
        void getIdentification(PeptideIdentification & id, MSSpectrum cid_spec, MSSpectrum etd_spec) nogil except + # wrap-doc:Performs and returns de novo identifications
