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

        AMSE_AdductInfo(AMSE_AdductInfo) nogil except + #wrap-ignore
        AMSE_AdductInfo(String & name, EmpiricalFormula & adduct, int charge, UInt mol_multiplier) nogil except +

        double getNeutralMass(double observed_mz) nogil except +
        double getMZ(double neutral_mass) nogil except +
        bool isCompatible(EmpiricalFormula db_entry) nogil except +
        int getCharge() nogil except +
        String getName() nogil except +
        # AMSE_AdductInfo parseAdductString(String & adduct) nogil except +

