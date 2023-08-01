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

        MSSim()      except + nogil 
        MSSim(MSSim &) except + nogil 

        void simulate(shared_ptr[SimRandomNumberGenerator] rnd_gen, SampleChannels peptides) except + nogil 
            # wrap-doc:
                #  General purpose function to simulate a mass spectrometry run
                #  
                #  
                #  :param rnd_gen: Random number generator which will be passed to the different classes
                #  :param peptides: List of peptides and abundances that will be simulated
        
        MSExperiment getExperiment() except + nogil  # wrap-doc:Returns the simulated experiment

        FeatureMap getSimulatedFeatures() except + nogil  # wrap-doc:Returns the simulated features

        ConsensusMap getChargeConsensus() except + nogil  # wrap-doc:Returns the charge consensus map of simulated features

        FeatureMap getContaminants() except + nogil  # wrap-doc:Returns the contaminants feature map of simulated features

        ConsensusMap getLabelingConsensus() except + nogil  # wrap-doc:Returns the labeling consensus map of simulated features

        MSExperiment getPeakMap() except + nogil  # wrap-doc:Returns the labeling consensus map of simulated features

        Param getParameters() except + nogil  # wrap-doc:Returns the default parameters for simulation including the labeling technique with name `labeling_name`

        void getIdentifications(libcpp_vector[ ProteinIdentification ] & proteins, libcpp_vector[ PeptideIdentification ] & peptides) except + nogil  # wrap-doc:Returns the simulated identifications (proteins and peptides)

        void getMS2Identifications(libcpp_vector[ ProteinIdentification ] & proteins, libcpp_vector[ PeptideIdentification ] & peptides) except + nogil  # wrap-doc:Returns the simulated MS2 identifications (proteins and peptides)

        void getFeatureIdentifications(libcpp_vector[ ProteinIdentification ] & proteins, libcpp_vector[ PeptideIdentification ] & peptides) except + nogil  # wrap-doc:Returns the simulated identifications (proteins and peptides) from feature annotations

