from Types cimport *
from String cimport *
from OpenSwathDataStructures cimport *
from ISpectrumAccess cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSCached.h>" namespace "OpenMS":

  # TODO missing functions
  cdef cppclass SpectrumAccessOpenMSCached(ISpectrumAccess):
        # wrap-inherits:
        #  ISpectrumAccess

        SpectrumAccessOpenMSCached() # wrap-pass-constructor

        SpectrumAccessOpenMSCached(String filename) nogil except +
        SpectrumAccessOpenMSCached(SpectrumAccessOpenMSCached q) nogil except + # wrap-ignore

