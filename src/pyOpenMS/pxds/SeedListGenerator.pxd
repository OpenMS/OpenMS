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
from DPosition cimport DPosition2


ctypedef libcpp_vector[DPosition2] SeedList

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/SeedListGenerator.h>" namespace "OpenMS":

    cdef cppclass SeedListGenerator:

        SeedListGenerator() nogil except +
        SeedListGenerator(SeedListGenerator &) nogil except + # compiler

        void generateSeedList(MSExperiment exp, libcpp_vector[DPosition2] & seeds) nogil except + # wrap-doc:Generate a seed list based on an MS experiment
        void generateSeedList(libcpp_vector[PeptideIdentification] & peptides, libcpp_vector[DPosition2] & seeds, bool use_peptide_mass) nogil except + # wrap-doc:Generates a seed list based on a list of peptide identifications
        # TODO map with UInt64
        void generateSeedList(ConsensusMap & consensus, Map[UInt64, libcpp_vector[DPosition2] ] & seeds) nogil except +  # wrap-ignore

        # TODO nested STL
        # void generateSeedLists(ConsensusMap & consensus, Map[ UInt64, libcpp_vector[ DPosition2] ] & seed_lists) nogil except + # wrap-doc:Generates seed lists based on a consensus map
        void convertSeedList(libcpp_vector[ DPosition2] & seeds, FeatureMap & features) nogil except + # wrap-doc:Converts a list of seed positions to a feature map (expected format for FeatureFinder)
        void convertSeedList(FeatureMap & features, libcpp_vector[ DPosition2] & seeds) nogil except + # wrap-doc:Converts a feature map with seed positions back to a simple list
