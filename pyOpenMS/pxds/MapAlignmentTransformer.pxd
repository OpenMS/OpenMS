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

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>" namespace "OpenMS":

    cdef cppclass MapAlignmentTransformer:

        MapAlignmentTransformer() nogil except +
        # /// Applies the <i>given</i> transformations to peak maps
        void transformPeakMaps(libcpp_vector[MSExperiment[Peak1D,ChromatogramPeak] ] & maps, libcpp_vector[TransformationDescription] & given_trafos)

        # /// Applies the <i>given</i> transformations to feature maps
        void transformFeatureMaps(libcpp_vector[FeatureMap[Feature] ] & maps, libcpp_vector[TransformationDescription] & given_trafos)

        # /// Applies the <i>given</i> transformations to consensus maps
        void transformConsensusMaps(libcpp_vector[ConsensusMap] & maps, libcpp_vector[TransformationDescription] & given_trafos)

        # NESTED std::vector TODO
        # /// Applies the <i>given</i> transformations to peptide identifications
        ## void transformPeptideIdentifications(libcpp_vector[libcpp_vector[PeptideIdentification] ] & maps, libcpp_vector[TransformationDescription] & given_trafos)

        # /// Applies the <i>given</i> transformations to a single peak map
        void transformSinglePeakMap(MSExperiment[Peak1D,ChromatogramPeak] & msexp, TransformationDescription & trafo)

        # /// Applies the <i>given</i> transformations to a single feature map
        void transformSingleFeatureMap(FeatureMap[Feature] & fmap, TransformationDescription & trafo)

        # /// Applies the <i>given</i> transformations to a single consensus map
        void transformSingleConsensusMap(ConsensusMap & cmap, TransformationDescription & trafo)

        # /// Applies the <i>given</i> transformations to a single peptide identification
        void transformSinglePeptideIdentification(libcpp_vector[PeptideIdentification] & pepids, TransformationDescription & trafo)

