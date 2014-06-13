from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from TransformationDescription cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from ConsensusMap cimport *
from FeatureMap cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>" namespace "OpenMS":
    
    cdef cppclass MapAlignmentAlgorithm(DefaultParamHandler,ProgressLogger) :
        # wrap-inherits:
        #  DefaultParamHandler
        #  ProgressLogger
        MapAlignmentAlgorithm() nogil except +
        MapAlignmentAlgorithm(MapAlignmentAlgorithm) nogil except + #wrap-ignore
        void alignPeakMaps(libcpp_vector[ MSExperiment[Peak1D, ChromatogramPeak] ] & , libcpp_vector[ TransformationDescription ] & ) nogil except +
        # TODO nested STL
        # void alignCompactFeatureMaps(libcpp_vector[ libcpp_vector[ Peak2D ] ] & , libcpp_vector[ TransformationDescription ] & ) nogil except +
        void alignFeatureMaps(libcpp_vector[ FeatureMap[Feature] ] & , libcpp_vector[ TransformationDescription ] & ) nogil except +
        void alignConsensusMaps(libcpp_vector[ ConsensusMap ] & , libcpp_vector[ TransformationDescription ] & ) nogil except +
        # void alignPeptideIdentifications(libcpp_vector[ libcpp_vector[ PeptideIdentification ] ] & , libcpp_vector[ TransformationDescription ] & ) nogil except +
        void setReference(Size reference_index, String & reference_file) nogil except +
        void fitModel(String & model_type, Param & params, libcpp_vector[ TransformationDescription ] & trafos) nogil except +
        void registerChildren() nogil except +

