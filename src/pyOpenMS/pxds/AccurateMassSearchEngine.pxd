from Types cimport *
from MassTrace cimport *
from Feature cimport *
from ConsensusFeature cimport *
from ConsensusMap cimport *
from FeatureMap cimport *
from MzTab cimport *
from DefaultParamHandler cimport *
from String cimport *
from ProgressLogger cimport *
from AccurateMassSearchResult cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>" namespace "OpenMS":
    
    cdef cppclass AccurateMassSearchEngine(DefaultParamHandler,ProgressLogger) :
        # wrap-inherits:
        #  DefaultParamHandler
        #  ProgressLogger
        AccurateMassSearchEngine() nogil except +
        AccurateMassSearchEngine(AccurateMassSearchEngine) nogil except + #wrap-ignore

        void queryByMass(double adduct_mass, Int adduct_charge , libcpp_vector[ AccurateMassSearchResult ] & ) nogil except +
        void queryByFeature(Feature feature, Size feature_index, libcpp_vector[ AccurateMassSearchResult ] & ) nogil except +
        void queryByConsensusFeature(ConsensusFeature cfeat, Size cf_index, Size number_of_maps, libcpp_vector[AccurateMassSearchResult]& results) nogil except +

        void run(FeatureMap & , MzTab & ) nogil except +
        void run(ConsensusMap&, MzTab&) nogil except +

        String getInternalIonMode() nogil except +

