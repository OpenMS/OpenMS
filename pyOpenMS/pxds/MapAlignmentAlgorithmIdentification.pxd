from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from Param cimport *
from Feature cimport *
from FeatureMap cimport *
from ConsensusMap cimport *
from TransformationDescription cimport *
from PeptideIdentification cimport *

from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>" namespace "OpenMS":

    cdef cppclass MapAlignmentAlgorithmIdentification:

        MapAlignmentAlgorithmIdentification() nogil except +
        void alignPeakMaps(libcpp_vector[MSExperiment[Peak1D,ChromatogramPeak]] & maps, libcpp_vector[TransformationDescription] & trafos) nogil except +
        void alignFeatureMaps(libcpp_vector[FeatureMap[Feature]] & features, libcpp_vector[TransformationDescription] & trafos) nogil except +
        void alignConsensusMaps(libcpp_vector[ConsensusMap] & features, libcpp_vector[TransformationDescription] & trafos) nogil except +
        # TODO nested STL
        void alignPeptideIdentifications(libcpp_vector[libcpp_vector[PeptideIdentification]] & ids, libcpp_vector[TransformationDescription] & trafos) nogil except + #wrap-ignore

        void setReference(Size ref, String model_type) nogil except +
        void fitModel(String model_type, Param & params, libcpp_vector[TransformationDescription] & trafos) nogil except +

