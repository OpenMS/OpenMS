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
        # wrap-doc:
                #   An implementation of the Spectrum Access interface using on-disk caching
                #   -----
                #   This class implements the OpenSWATH Spectrum Access interface
                #   (ISpectrumAccess) using the CachedmzML class which is able to read and
                #   write a cached mzML file

        SpectrumAccessOpenMSCached(SpectrumAccessOpenMSCached &) nogil except +

