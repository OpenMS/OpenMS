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
        # wrap-doc:
        #  FLASHDeconv algorithm: ultrafast mass deconvolution algorithm for top down mass spectrometry dataset.
        #  From MSSpectrum, this class outputs DeconvolvedSpectrum.
        #  Deconvolution takes three steps:
        #    i) decharging and select candidate masses - speed up via binning
        #    ii) collecting isotopes from the candidate masses and deisotoping - peak groups are defined here
        #    iii) scoring and filter out low scoring masses (i.e., peak groups)

        # Constructors
        FLASHDeconvAlgorithm() except + nogil
        FLASHDeconvAlgorithm(FLASHDeconvAlgorithm &) except + nogil

        # Main run method
        void run(MSExperiment & input_map, libcpp_vector[DeconvolvedSpectrum] & deconvolved_spectra, libcpp_vector[MassFeature_FDHS] & deconvolved_features) except + nogil
        # wrap-doc:
        #  Run FLASHDeconv algorithm for input_map and store deconvolved_spectra and deconvolved_features.
        #  :param input_map: The input MSExperiment containing spectra to deconvolve
        #  :param deconvolved_spectra: Output vector to store deconvolved spectra
        #  :param deconvolved_features: Output vector to store mass features

        # Averagine access
        PrecalAveragine& getAveragine() except + nogil
        # wrap-doc:Get calculated averagine. Call after run() is called.

        PrecalAveragine& getDecoyAveragine() except + nogil
        # wrap-doc:Get calculated decoy averagine. Call after run() is called.

        # Tolerances
        libcpp_vector[double] getTolerances() except + nogil
        # wrap-doc:Get mass tolerances per MS level.

        # Noise decoy weight
        double getNoiseDecoyWeight() except + nogil
        # wrap-doc:Get noise decoy weight determined during q-value calculation.

        # Static method
        # Note: Static methods need special handling in Cython


cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>" namespace "OpenMS::FLASHDeconvAlgorithm":

    int getScanNumber(MSExperiment & exp, Size index) except + nogil  # wrap-attach:FLASHDeconvAlgorithm
    # wrap-doc:Get scan number of the spectrum at index in exp
