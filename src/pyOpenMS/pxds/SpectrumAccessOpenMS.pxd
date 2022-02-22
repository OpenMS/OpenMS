from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from OpenSwathDataStructures cimport *
from ISpectrumAccess cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>" namespace "OpenMS":

  cdef cppclass SpectrumAccessOpenMS(ISpectrumAccess):
        # wrap-inherits:
        #  ISpectrumAccess
        # wrap-doc:
        #   An implementation of the OpenSWATH Spectrum Access interface using OpenMS

        SpectrumAccessOpenMS() # wrap-pass-constructor

        SpectrumAccessOpenMS(SpectrumAccessOpenMS &) nogil except + 
        SpectrumAccessOpenMS(shared_ptr[ MSExperiment ] & ms_experiment) nogil except +
