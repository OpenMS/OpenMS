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
        OfflinePrecursorIonSelection(OfflinePrecursorIonSelection &) nogil except + # compiler

        void makePrecursorSelectionForKnownLCMSMap(FeatureMap & features,
                                                   MSExperiment & experiment,
                                                   MSExperiment & ms2,
                                                   libcpp_set[ int ] & charges_set, bool
                                                   feature_based) nogil except +
            # wrap-doc:
                #   Makes the precursor selection for a given feature map, either feature or scan based
                #   -----
                #   :param features: Input feature map
                #   :param experiment: Input raw data
                #   :param ms2: Precursors are added as empty MS2 spectra to this MSExperiment
                #   :param charges_set: Allowed charge states
                #   :param feature_based: If true the selection is feature based, if false it is scan based and the highest signals in each spectrum are chosen

        # TODO nested STL
        void getMassRanges(FeatureMap & features, 
                           MSExperiment & experiment,
                           libcpp_vector[ libcpp_vector[ libcpp_pair[ Size, Size ] ] ] & indices) nogil except + # wrap-ignore

        void createProteinSequenceBasedLPInclusionList(String include_, String rt_model_file, String pt_model_file, FeatureMap & precursors) nogil except +
        void setLPSolver(SOLVER solver) nogil except +
        SOLVER getLPSolver() nogil except +

