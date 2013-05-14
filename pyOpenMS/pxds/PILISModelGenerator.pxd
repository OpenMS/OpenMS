from Types cimport *
from DefaultParamHandler cimport *
from HiddenMarkovModel cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/PILISModelGenerator.h>" namespace "OpenMS":
    
    cdef cppclass PILISModelGenerator(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        PILISModelGenerator() nogil except +
        PILISModelGenerator(PILISModelGenerator) nogil except +
        void getModel(HiddenMarkovModel & model) nogil except +

