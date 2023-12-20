# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/ParentPeakMower.h>" namespace "OpenMS":

    cdef cppclass ParentPeakMower(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler
        # wrap-doc:
        #  ParentPeakMower gets rid of high peaks that could stem from unfragmented precursor ions
        
        ParentPeakMower() except + nogil 
        ParentPeakMower(ParentPeakMower &) except + nogil  

        void filterSpectrum(MSSpectrum & spec) except + nogil 
        void filterPeakSpectrum(MSSpectrum & spec) except + nogil 
        void filterPeakMap(MSExperiment & exp) except + nogil 
