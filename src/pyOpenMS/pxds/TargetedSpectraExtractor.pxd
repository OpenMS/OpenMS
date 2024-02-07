from TransitionTSVFile cimport *
from TargetedExperiment cimport *
from DataValue cimport *
from DefaultParamHandler cimport *
from MSExperiment cimport *
from FeatureMap cimport *
from String cimport *
from Types cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.h>" namespace "OpenMS":

    cdef cppclass TargetedSpectraExtractor(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler

        TargetedSpectraExtractor() except + nogil 
        TargetedSpectraExtractor(TargetedSpectraExtractor &) except + nogil  # compiler

        void getDefaultParameters(Param&) except + nogil 

        void annotateSpectra(libcpp_vector[ MSSpectrum ]&, TargetedExperiment&, libcpp_vector[ MSSpectrum ]&, FeatureMap&) except + nogil 
        void annotateSpectra(libcpp_vector[ MSSpectrum ]&, TargetedExperiment&, libcpp_vector[ MSSpectrum ]&) except + nogil 
        void annotateSpectra(libcpp_vector[ MSSpectrum ]&, FeatureMap&, FeatureMap&, libcpp_vector[ MSSpectrum ]&) except + nogil 

        void searchSpectrum(FeatureMap&, FeatureMap&, bool) except + nogil 

        void pickSpectrum(MSSpectrum&, MSSpectrum&) except + nogil 

        void scoreSpectra(libcpp_vector[ MSSpectrum ]&, libcpp_vector[ MSSpectrum ]&, FeatureMap&, libcpp_vector[ MSSpectrum ]&) except + nogil 
        void scoreSpectra(libcpp_vector[ MSSpectrum ]&, libcpp_vector[ MSSpectrum ]&, libcpp_vector[ MSSpectrum ]&) except + nogil 

        void selectSpectra(libcpp_vector[ MSSpectrum ]&, FeatureMap&, libcpp_vector[ MSSpectrum ]&, FeatureMap&) except + nogil 
        void selectSpectra(libcpp_vector[ MSSpectrum ]&, libcpp_vector[ MSSpectrum ]&) except + nogil 

        void extractSpectra(MSExperiment&, TargetedExperiment&, libcpp_vector[ MSSpectrum ]&, FeatureMap&, bool) except + nogil 
        void extractSpectra(MSExperiment&, TargetedExperiment&, libcpp_vector[ MSSpectrum ]&) except + nogil 
        void extractSpectra(MSExperiment&, FeatureMap&, libcpp_vector[ MSSpectrum ]&) except + nogil 

        void constructTransitionsList(FeatureMap&, FeatureMap&, TargetedExperiment&) except + nogil 

        void storeSpectraMSP(const String&, MSExperiment&) except + nogil 

        void mergeFeatures(FeatureMap&, FeatureMap&) except + nogil 

        # void matchSpectrum(MSSpectrum& input_spectrum, TSE_Comparator& cmp, libcpp_vector[ TSE_Match ]& matches) except + nogil 


cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.h>" namespace "OpenMS::TargetedSpectraExtractor":

    cdef cppclass TSE_Match "OpenMS::TargetedSpectraExtractor::Match":

        TSE_Match() except + nogil 
        TSE_Match(TSE_Match &) except + nogil 
        TSE_Match(MSSpectrum& spectrum, double score) except + nogil 

        MSSpectrum spectrum
        double score

    # cdef cppclass TSE_Comparator "OpenMS::TargetedSpectraExtractor::Comparator":

        # TSE_Comparator() except + nogil 
        # TSE_Comparator(TSE_Comparator &) except + nogil 

        # void generateScores(MSSpectrum& spec, libcpp_vector[libcpp_pair[Size,double]]& scores, double min_score) except + nogil 
        # void init(libcpp_vector[MSSpectrum]& library, libcpp_map[String,DataValue]& options) except + nogil 
        # libcpp_vector[MSSpectrum]& getLibrary() except + nogil 

    # cdef cppclass TSE_BinnedSpectrumComparator "OpenMS::TargetedSpectraExtractor::BinnedSpectrumComparator" (TSE_Comparator):
        # wrap-inherits:
        #  TSE_Comparator

        # TSE_BinnedSpectrumComparator() except + nogil 
        # TSE_BinnedSpectrumComparator(TSE_BinnedSpectrumComparator &) except + nogil 

        # void generateScores(MSSpectrum& spec, libcpp_vector[libcpp_pair[Size,double]]& scores, double min_score) except + nogil 
        # void init(libcpp_vector[MSSpectrum]& library, libcpp_map[String,DataValue]& options) except + nogil 
