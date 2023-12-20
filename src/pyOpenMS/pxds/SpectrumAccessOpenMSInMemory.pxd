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

        # 'ISpectrumAccess' is not a type identifier
        SpectrumAccessOpenMSInMemory(ISpectrumAccess &) except + nogil  # wrap-ignore
        SpectrumAccessOpenMSInMemory(SpectrumAccessOpenMS &) except + nogil 
        SpectrumAccessOpenMSInMemory(SpectrumAccessOpenMSCached &) except + nogil 
        SpectrumAccessOpenMSInMemory(SpectrumAccessOpenMSInMemory &) except + nogil 
        SpectrumAccessOpenMSInMemory(SpectrumAccessQuadMZTransforming &) except + nogil 
