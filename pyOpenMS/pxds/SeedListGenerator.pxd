from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from MSSpectrum cimport *
from Map cimport *
from ConsensusMap cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from PeptideIdentification cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *


ctypedef libcpp_vector[ DPosition2] SeedList

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/SeedListGenerator.h>" namespace "OpenMS":

    cdef cppclass SeedListGenerator:

        SeedListGenerator()                    nogil except +
        SeedListGenerator(SeedListGenerator &) nogil except +

        void generateSeedList(MSExperiment[Peak1D, ChromatogramPeak] exp, libcpp_vector[DPosition2] & seeds) nogil except +
        void generateSeedList(libcpp_vector[PeptideIdentification] & peptides, libcpp_vector[DPosition2] & seeds, bool use_peptide_mass) nogil except +
        # TODO map with UInt64
        void generateSeedList(ConsensusMap & consensus, Map[unsigned long, libcpp_vector[DPosition2] ] & seeds) nogil except +  # wrap-ignore
        # TODO nested STL
        # void generateSeedLists(ConsensusMap & consensus, Map[ UInt64, libcpp_vector[ DPosition2] ] & seed_lists)
        void convertSeedList(libcpp_vector[ DPosition2] & seeds, FeatureMap[Feature] & features)
        void convertSeedList(FeatureMap[Feature] & features, libcpp_vector[ DPosition2] & seeds)
