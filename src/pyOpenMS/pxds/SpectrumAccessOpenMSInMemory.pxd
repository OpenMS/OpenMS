from Types cimport *
from OpenSwathDataStructures cimport *
from MSSpectrum cimport *
from MSChromatogram cimport *
from MSExperiment cimport *
from ISpectrumAccess cimport *

from SpectrumAccessOpenMS cimport *
from SpectrumAccessOpenMSCached cimport *
from SpectrumAccessQuadMZTransforming cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>" namespace "OpenMS":
    
    cdef cppclass SpectrumAccessOpenMSInMemory(ISpectrumAccess) :
        # wrap-inherits:
        #  ISpectrumAccess

        SpectrumAccessOpenMSInMemory() # wrap-pass-constructor

        # SpectrumAccessOpenMSInMemory(ISpectrumAccess & origin) nogil except +
        SpectrumAccessOpenMSInMemory(SpectrumAccessOpenMS) nogil except +
        SpectrumAccessOpenMSInMemory(SpectrumAccessOpenMSCached) nogil except +
        SpectrumAccessOpenMSInMemory(SpectrumAccessOpenMSInMemory) nogil except +
        SpectrumAccessOpenMSInMemory(SpectrumAccessQuadMZTransforming) nogil except +

