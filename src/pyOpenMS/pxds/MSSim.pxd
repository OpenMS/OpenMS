from libcpp.map cimport map as libcpp_map
from libcpp cimport bool
from smart_ptr cimport shared_ptr
from Types cimport *
from ProgressLogger cimport *
from DefaultParamHandler cimport *
from Param cimport *
from SimTypes cimport *
from ConsensusMap cimport *
from MSExperiment cimport *
from FeatureMap cimport *

# SimTypes::MSSimExperiment = MSExperiment
# SimTypes::FeatureMapSim = FeatureMap
cdef extern from "<OpenMS/SIMULATION/MSSim.h>" namespace "OpenMS":

    cdef cppclass MSSim:

        MSSim()      nogil except +
        MSSim(MSSim) nogil except + # wrap-ignore

        # General purpose function to simulate a mass spectrometry run
        #
        #@param rnd_gen random number generator which will be passed to the different classes
        #@param peptides List of peptides and abundances that will be simulated
        void simulate(shared_ptr[SimRandomNumberGenerator] rnd_gen, SampleChannels peptides) nogil except +

        # Access the simulated experiment
        MSExperiment getExperiment() nogil except +

        # Access the simulated features
        FeatureMap getSimulatedFeatures() nogil except +

        # Access the charge consensus map of simulated features
        ConsensusMap getChargeConsensus() nogil except +

        # Access the contaminants feature map of simulated features
        FeatureMap getContaminants() nogil except +

        # Access the labeling consensus map of simulated features
        ConsensusMap getLabelingConsensus() nogil except +

        # Access the labeling consensus map of simulated features
        MSExperiment getPeakMap() nogil except +

        # Returns the default parameters for simulation including the labeling technique with name @p labeling_name
        Param getParameters() nogil except +

        void getIdentifications(libcpp_vector[ ProteinIdentification ] & proteins, libcpp_vector[ PeptideIdentification ] & peptides) nogil except +

        void getMS2Identifications(libcpp_vector[ ProteinIdentification ] & proteins, libcpp_vector[ PeptideIdentification ] & peptides) nogil except +

        void getFeatureIdentifications(libcpp_vector[ ProteinIdentification ] & proteins, libcpp_vector[ PeptideIdentification ] & peptides) nogil except +

