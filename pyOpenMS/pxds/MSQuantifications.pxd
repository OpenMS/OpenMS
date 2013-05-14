from libcpp cimport vector as libcpp_vector
from libcpp cimport map as libcpp_map
from libcpp cimport pair as libcpp_pair
from ConsensusMap cimport *
from String cimport *
from StringList cimport *
from ExperimentalSettings cimport *
from DataProcessing cimport *
from MetaInfo cimport *
from ConsensusMap cimport *
from CVTermList cimport *
from FeatureMap cimport *
from Feature cimport *

cdef extern from "<OpenMS/METADATA/MSQuantifications.h>" namespace "OpenMS":

    cdef cppclass MSQuantifications:
        MSQuantifications()  nogil
        MSQuantifications(MSQuantifications) nogil

        bool operator==(MSQuantifications &) nogil
        bool operator!=(MSQuantifications &) nogil

        # not yet implemented
        # void load(String filename, bool trim_lines, Int first_n) nogil except +

        libcpp_vector[DataProcessing] getDataProcessingList() nogil except +
        libcpp_vector[Assay] getAssays() nogil except +
        # TODO STL map with wrapped key
        # libcpp_map[String, Ratio] getRatios() nogil except +
        libcpp_vector[ConsensusMap] getConsensusMaps() nogil except +
        void setConsensusMaps(libcpp_vector[ConsensusMap]) nogil except +
        libcpp_vector[FeatureMap[Feature] ] getFeatureMaps() nogil except +
        AnalysisSummary getAnalysisSummary() nogil except +
        void setDataProcessingList(libcpp_vector[DataProcessing] dpl) nogil except +
        void setAnalysisSummaryQuantType(QUANT_TYPES r) nogil except +
        void addConsensusMap(ConsensusMap m) nogil except +
        void assignUIDs() nogil except +
        # TODO nested STL
        # void registerExperiment(MSExperiment[Peak1D] exp, libcpp_vector[libcpp_vector[libcpp_pair[String, DoubleReal] ] ] labels) nogil except +

cdef extern from "<OpenMS/METADATA/MSQuantifications.h>" namespace "OpenMS::MSQuantifications":
    # derived from processing applied
    cdef enum QUANT_TYPES:
        # wrap-attach:
        #     MSQuantifications
        MS1LABEL = 0,
        MS2LABEL,
        LABELFREE,
        SIZE_OF_QUANT_TYPES

    cdef cppclass AnalysisSummary:
        # wrap-attach:
        #     MSQuantifications
        AnalysisSummary() 
        AnalysisSummary(AnalysisSummary)

        MetaInfo user_params_
        CVTermList cv_params_
        QUANT_TYPES quant_type_

    cdef cppclass Assay:
        # wrap-attach:
        #     MSQuantifications
        Assay()
        Assay(Assay)

        String uid_ 
        # libcpp_vector[libcpp_pair[String, double] ] mods_
        libcpp_vector[ExperimentalSettings] raw_files_
        libcpp_map[size_t, FeatureMap[Feature] ] feature_maps_ 
        # iTRAQ needs no FeatureMaps so ExperimentalSettings are not directly mapped to FeatureMaps

