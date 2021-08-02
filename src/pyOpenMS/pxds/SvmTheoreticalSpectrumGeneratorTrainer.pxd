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
        #    DefaultParamHandler

        SvmTheoreticalSpectrumGeneratorTrainer() nogil except +
        SvmTheoreticalSpectrumGeneratorTrainer(SvmTheoreticalSpectrumGeneratorTrainer &) nogil except +

        void trainModel(MSExperiment & spectra, libcpp_vector[AASequence] & annotations, String filename, int precursor_charge) nogil except +
        void normalizeIntensity(MSSpectrum & S) nogil except +
