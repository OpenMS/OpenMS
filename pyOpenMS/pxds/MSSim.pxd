from libcpp.map cimport map as libcpp_map
from libcpp cimport bool
from Types cimport *
from ProgressLogger cimport *
from DefaultParamHandler cimport *
from Param cimport *
from SimTypes cimport *
from ConsensusMap cimport *
from MSExperiment cimport *
from FeatureMap cimport *

# MSSimExperiment = MSExperiment
# FeatureMapSim = FeatureMap[Feature]
cdef extern from "<OpenMS/SIMULATION/MSSim.h>" namespace "OpenMS":

    cdef cppclass MSSim:

        MSSim()      nogil except +
        MSSim(MSSim) nogil except + # wrap-ignore

        # General purpose function to simulate a mass spectrometry run
        #
        #@param rnd_gen GSL random number generator which will be passed to the different classes
        #@param peptides List of peptides and abundances that will be simulated
        void simulate(SimRandomNumberGenerator rnd_gen, SampleChannels peptides) nogil except +

        # Access the simulated experiment
        MSExperiment[Peak1D, ChromatogramPeak] getExperiment() nogil except +

        # Access the simulated features
        FeatureMap[Feature] getSimulatedFeatures() nogil except +

        # Access the charge consensus map of simulated features
        ConsensusMap getChargeConsensus() nogil except +

        # Access the contaminants feature map of simulated features
        FeatureMap[Feature] getContaminants() nogil except +

        # Access the labeling consensus map of simulated features
        ConsensusMap getLabelingConsensus() nogil except +

        # Access the labeling consensus map of simulated features
        MSExperiment[Peak1D, ChromatogramPeak] getPeakMap() nogil except +

        # Returns the default parameters for simulation including the labeling technique with name @p labeling_name
        Param getParameters() nogil except +
