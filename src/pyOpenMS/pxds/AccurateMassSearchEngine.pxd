from Types cimport *
from MassTrace cimport *
from MzTabM cimport *
from Feature cimport *
from ConsensusFeature cimport *
from ConsensusMap cimport *
from FeatureMap cimport *
from MzTab cimport *
from DefaultParamHandler cimport *
from String cimport *
from ProgressLogger cimport *
from AccurateMassSearchResult cimport *
from EmpiricalFormula cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>" namespace "OpenMS":

    cdef cppclass AccurateMassSearchEngine(DefaultParamHandler,ProgressLogger) :
        # wrap-inherits:
        #  DefaultParamHandler
        #  ProgressLogger
        AccurateMassSearchEngine() nogil except +
        AccurateMassSearchEngine(AccurateMassSearchEngine &) nogil except +

        void queryByMZ(double observed_mz, Int observed_charge, String ion_mode,
                         libcpp_vector[ AccurateMassSearchResult ] &, EmpiricalFormula & observed_adduct) nogil except +

        void queryByFeature(Feature feature, Size feature_index, String ion_mode,
                            libcpp_vector[ AccurateMassSearchResult ] & ) nogil except +

        void queryByConsensusFeature(ConsensusFeature cfeat, Size cf_index, Size number_of_maps, String ion_mode,
                                     libcpp_vector[AccurateMassSearchResult]& results) nogil except +

        void run(FeatureMap&, MzTab&) nogil except +
        void run(FeatureMap&, MzTabM&) nogil except +
        void run(ConsensusMap&, MzTab&) nogil except +

        void init() nogil except +


## cdef extern from "<OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>" namespace "OpenMS":
##     
##     cdef cppclass AdductInfo "OpenMS::AdductInfo":
## 
##         #private
##         AdductInfo() nogil except + # wrap-ignore
##         #private
##         AdductInfo(AdductInfo) nogil except + # wrap-ignore 
##         AdductInfo(String name, EmpiricalFormula & adduct, int charge, UInt mol_multiplier) nogil except +
##         double getNeutralMass(double observed_mz) nogil except +
##         double getMZ(double neutral_mass) nogil except +
##         bool isCompatible(EmpiricalFormula db_entry) nogil except +
##         int getCharge() nogil except +
##         String getName() nogil except +
##         UInt getMolMultiplier() nogil except +
