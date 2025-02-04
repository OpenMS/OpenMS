from Types cimport *
from MassTrace cimport *
from Feature cimport *
from ConsensusFeature cimport *
from FeatureMap cimport *
from ConsensusMap cimport *
from MzTab cimport *
from DefaultParamHandler cimport *
from String cimport *
from EmpiricalFormula cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>" namespace "OpenMS":
    
    cdef cppclass AMSE_AdductInfo "OpenMS::AdductInfo":

        # private
        AMSE_AdductInfo() except + nogil  #wrap-ignore
        # private
        AMSE_AdductInfo(AMSE_AdductInfo) except + nogil  #wrap-ignore

        AMSE_AdductInfo(const String & name, EmpiricalFormula & adduct, int charge, UInt mol_multiplier) except + nogil 

        double getNeutralMass(double observed_mz) except + nogil 
        double getMZ(double neutral_mass) except + nogil 
        bool isCompatible(EmpiricalFormula db_entry) except + nogil 
        int getCharge() except + nogil 
        String getName() except + nogil 

        # AMSE_AdductInfo parseAdductString(const String & adduct) except + nogil 

