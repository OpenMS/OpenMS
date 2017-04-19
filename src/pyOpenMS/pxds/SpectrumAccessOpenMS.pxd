from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from OpenSwathDataStructures cimport *
from ISpectrumAccess cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>" namespace "OpenMS":

  cdef cppclass SpectrumAccessOpenMS(ISpectrumAccess):
        # wrap-inherits:
        #  ISpectrumAccess

        SpectrumAccessOpenMS() # wrap-pass-constructor

        SpectrumAccessOpenMS(SpectrumAccessOpenMS)
        SpectrumAccessOpenMS(shared_ptr[ MSExperiment ] & ms_experiment) nogil except +

