from libcpp.vector cimport vector as libcpp_vector
from PeptideIdentification cimport *
from ProteinIdentification cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>" namespace "OpenMS":

    cdef cppclass FalseDiscoveryRate(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        FalseDiscoveryRate()                  nogil except +
        FalseDiscoveryRate(FalseDiscoveryRate)   nogil except + #wrap-ignore

        void apply(libcpp_vector[PeptideIdentification] & forward_ids, libcpp_vector[PeptideIdentification] & reverse_ids) nogil except + 
        void apply(libcpp_vector[PeptideIdentification] & id) nogil except + 
        void apply(libcpp_vector[ProteinIdentification] & forward_ids, libcpp_vector[ProteinIdentification] & reverse_ids) nogil except + 
        void apply(libcpp_vector[ProteinIdentification] & id) nogil except + 

