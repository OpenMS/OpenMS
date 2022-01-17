# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from String cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>" namespace "OpenMS":

    cdef cppclass SpectraMerger(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        SpectraMerger() nogil except + # wrap-doc:Merges blocks of MS or MS2 spectra
        SpectraMerger(SpectraMerger &) nogil except +

        void mergeSpectraBlockWise(MSExperiment & exp) nogil except +
        void mergeSpectraPrecursors(MSExperiment & exp) nogil except + # wrap-doc:Merges spectra with similar precursors (must have MS2 level)
        void average(MSExperiment & exp, String average_type) nogil except +
        # wrap-doc:
                #   Average over neighbouring spectra
                #   -----
                #   :param exp: Experimental data to be averaged
                #   :param average_type: Averaging type to be used ("gaussian" or "tophat")
