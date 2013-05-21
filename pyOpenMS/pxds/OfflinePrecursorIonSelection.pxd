from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
from libcpp.set cimport set as libcpp_set
from libcpp.vector cimport vector as libcpp_vector
from FeatureMap cimport *
from MSExperiment cimport *
from IDMapper cimport *
from FeatureXMLFile cimport *
# from PSLPFormulation cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from LPWrapper cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/OfflinePrecursorIonSelection.h>" namespace "OpenMS":
    
    cdef cppclass OfflinePrecursorIonSelection(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        OfflinePrecursorIonSelection() nogil except +
        OfflinePrecursorIonSelection(OfflinePrecursorIonSelection) nogil except + #wrap-ignore
        void makePrecursorSelectionForKnownLCMSMap(FeatureMap[Feature] & features, MSExperiment[ Peak1D, ChromatogramPeak ] & experiment, MSExperiment[ Peak1D, ChromatogramPeak ] & ms2, libcpp_set[ int ] & charges_set, bool feature_based) nogil except +
        # TODO nested STL
        void getMassRanges(FeatureMap[Feature] & features, 
                           MSExperiment[ Peak1D, ChromatogramPeak ] & experiment,
                           libcpp_vector[ libcpp_vector[ libcpp_pair[ Size, Size ] ] ] & indices) nogil except + # wrap-ignore
        void createProteinSequenceBasedLPInclusionList(String include_, String rt_model_file, String pt_model_file, FeatureMap[Feature] & precursors) nogil except +
        void setLPSolver(SOLVER solver) nogil except +
        SOLVER getLPSolver() nogil except +

