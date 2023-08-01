from Types cimport *

from MSExperiment cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from AASequence cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGeneratorTrainer.h>" namespace "OpenMS":

    cdef cppclass SvmTheoreticalSpectrumGeneratorTrainer(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler

        SvmTheoreticalSpectrumGeneratorTrainer() except + nogil 
        SvmTheoreticalSpectrumGeneratorTrainer(SvmTheoreticalSpectrumGeneratorTrainer &) except + nogil 

        void trainModel(MSExperiment & spectra, libcpp_vector[AASequence] & annotations, String filename, int precursor_charge) except + nogil 
        void normalizeIntensity(MSSpectrum & S) except + nogil 
