from libcpp.vector cimport vector as libcpp_vector
from MSExperiment cimport *
from MSSpectrum cimport *
from Matrix cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from DeconvolvedSpectrum cimport *
from FLASHHelperClasses cimport *



cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>" namespace "OpenMS":

    cdef cppclass FLASHDeconvAlgorithm(DefaultParamHandler,ProgressLogger):
        # wrap-inherits:
        #   DefaultParamHandler
        #   ProgressLogger

        # default constructor
        FLASHDeconvAlgorithm() except + nogil
        # copy constructor
        FLASHDeconvAlgorithm(FLASHDeconvAlgorithm &) except + nogil
        # move constructor
        #FLASHDeconvAlgorithm(FLASHDeconvAlgorithm && other) except + nogil
        
        PrecalAveragine& getAveragine() except + nogil

        double getNoiseDecoyWeight() except + nogil

        
        
        void run(MSExperiment & input_map,libcpp_vector[DeconvolvedSpectrum] & Dspectrum, libcpp_vector[MassFeature_FDHS] & massfeature) except + nogil
