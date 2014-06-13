from Types cimport *
from DefaultParamHandler cimport *
from PeptideIdentification cimport *
from TextFile cimport *
from PeptideAndProteinQuant cimport *
from ProteinResolver cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/QuantitativeExperimentalDesign.h>" namespace "OpenMS":
    
    cdef cppclass QuantitativeExperimentalDesign(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        QuantitativeExperimentalDesign() nogil except +
        QuantitativeExperimentalDesign(QuantitativeExperimentalDesign) nogil except + #wrap-ignore
        void applyDesign2Resolver(ProteinResolver & resolver, TextFile & file_, StringList & fileNames) nogil except +
        void applyDesign2Quantifier(PeptideAndProteinQuant & quantifier, TextFile & file_, StringList & fileNames) nogil except +

