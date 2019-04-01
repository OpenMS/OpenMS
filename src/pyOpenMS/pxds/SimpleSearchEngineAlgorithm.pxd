from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from ProteinIdentification cimport *
from PeptideIdentification cimport *
from String cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/SimpleSearchEngineAlgorithm.h>" namespace "OpenMS":
 
    cdef cppclass SimpleSearchEngineAlgorithm(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        SimpleSearchEngineAlgorithm()   nogil except +
        SimpleSearchEngineAlgorithm(SimpleSearchEngineAlgorithm) nogil except + # wrap-ignore

        void search(const String& in_mzML, 
          const String& in_db, 
          libcpp_vector[ ProteinIdentification ] & prot_ids,
          libcpp_vector[ PeptideIdentification ] & pep_ids) nogil except +

