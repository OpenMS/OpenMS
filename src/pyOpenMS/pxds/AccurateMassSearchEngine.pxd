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
        AccurateMassSearchEngine() except + nogil 
        AccurateMassSearchEngine(AccurateMassSearchEngine &) except + nogil 

        void queryByMZ(double observed_mz, Int observed_charge, String ion_mode,
                         libcpp_vector[ AccurateMassSearchResult ] &, EmpiricalFormula & observed_adduct) except + nogil 

        void queryByFeature(Feature feature, Size feature_index, String ion_mode,
                            libcpp_vector[ AccurateMassSearchResult ] & ) except + nogil 

        void queryByConsensusFeature(ConsensusFeature cfeat, Size cf_index, Size number_of_maps, String ion_mode,
                                     libcpp_vector[AccurateMassSearchResult]& results) except + nogil 

        void run(FeatureMap&, MzTab&) except + nogil 
        void run(FeatureMap&, MzTabM&) except + nogil 
        void run(ConsensusMap&, MzTab&) except + nogil 

        void init() except + nogil 


## cdef extern from "<OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>" namespace "OpenMS":
##     
##     cdef cppclass AdductInfo "OpenMS::AdductInfo":
## 
##         #private
##         AdductInfo() except + nogil  # wrap-ignore
##         #private
##         AdductInfo(AdductInfo) except + nogil  # wrap-ignore 
##         AdductInfo(String name, EmpiricalFormula & adduct, int charge, UInt mol_multiplier) except + nogil 
##         double getNeutralMass(double observed_mz) except + nogil 
##         double getMZ(double neutral_mass) except + nogil 
##         bool isCompatible(EmpiricalFormula db_entry) except + nogil 
##         int getCharge() except + nogil 
##         String getName() except + nogil 
##         UInt getMolMultiplier() except + nogil 
