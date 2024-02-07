from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
from ConsensusMap cimport *
from String cimport *
from Types cimport *
from StringList cimport *
from ExperimentalSettings cimport *
from DataProcessing cimport *
from MetaInfo cimport *
from ConsensusMap cimport *
from CVTermList cimport *
from FeatureMap cimport *
from Feature cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/METADATA/MSQuantifications.h>" namespace "OpenMS":

    cdef cppclass MSQuantifications(ExperimentalSettings):
        # wrap-inherits:
        #  ExperimentalSettings
        MSQuantifications() except + nogil 
        MSQuantifications(MSQuantifications &) except + nogil 

        # Detailed constructor from a FeatureMap
        MSQuantifications(FeatureMap fm,
                          ExperimentalSettings& es,
                          libcpp_vector[DataProcessing]& dps) except + nogil 
        # std::vector<std::vector<std::pair<String, double> > > labels = (std::vector<std::vector<std::pair<String, double> > >()));
        
        bool operator==(MSQuantifications &) except + nogil 
        bool operator!=(MSQuantifications &) except + nogil 

        # TODO - not implemented in OpenMS, remove from API
        # void load(String filename, bool trim_lines, Int first_n) except + nogil 

        libcpp_vector[DataProcessing] getDataProcessingList() except + nogil 
        libcpp_vector[Assay] getAssays() except + nogil 

        # TODO - not implemented in OpenMS, remove from API
        # libcpp_map[String, Ratio] getRatios() except + nogil  # wrap-ignore

        libcpp_vector[ConsensusMap] getConsensusMaps() except + nogil 
        void setConsensusMaps(libcpp_vector[ConsensusMap]) except + nogil 
        libcpp_vector[FeatureMap ] getFeatureMaps() except + nogil 
        AnalysisSummary getAnalysisSummary() except + nogil 
        void setDataProcessingList(libcpp_vector[DataProcessing] dpl) except + nogil 
        void setAnalysisSummaryQuantType(QUANT_TYPES r) except + nogil 
        void addConsensusMap(ConsensusMap m) except + nogil 
        void assignUIDs() except + nogil 
        void registerExperiment(MSExperiment exp, 
                                libcpp_vector[ libcpp_vector[ libcpp_pair[
                                  String, double] ] ] labels) except + nogil  # wrap-ignore

cdef extern from "<OpenMS/METADATA/MSQuantifications.h>" namespace "OpenMS::MSQuantifications":
    # derived from processing applied
    cdef enum QUANT_TYPES:
        # wrap-attach:
        #    MSQuantifications
        MS1LABEL = 0,
        MS2LABEL,
        LABELFREE,
        SIZE_OF_QUANT_TYPES

    cdef cppclass AnalysisSummary:
        AnalysisSummary()  except + nogil 
        AnalysisSummary(AnalysisSummary &) except + nogil 

        MetaInfo user_params_
        CVTermList cv_params_
        QUANT_TYPES quant_type_

    cdef cppclass Assay:
        Assay() except + nogil 
        Assay(Assay &) except + nogil 

        String uid_ 
        libcpp_vector[libcpp_pair[String, double] ] mods_ # wrap-ignore
        libcpp_vector[ExperimentalSettings] raw_files_
        libcpp_map[size_t, FeatureMap ] feature_maps_ 
        # iTRAQ needs no FeatureMaps so ExperimentalSettings are not directly mapped to FeatureMaps

