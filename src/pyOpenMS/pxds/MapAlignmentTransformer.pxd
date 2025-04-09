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
        # wrap-doc:
            #  This class collects functions for applying retention time transformations to data structures

        MapAlignmentTransformer() except + nogil  
        MapAlignmentTransformer(MapAlignmentTransformer &) except + nogil 

        void transformRetentionTimes(MSExperiment&, TransformationDescription&, bool) except + nogil  # wrap-doc:Applies the given transformation to a peak map

        void transformRetentionTimes(FeatureMap&, TransformationDescription&, bool) except + nogil  # wrap-doc:Applies the given transformation to a feature map

        void transformRetentionTimes(ConsensusMap&, TransformationDescription&, bool) except + nogil  # wrap-doc:Applies the given transformation to a consensus map

        void transformRetentionTimes(libcpp_vector[PeptideIdentification]&, TransformationDescription&, bool) except + nogil  # wrap-doc:Applies the given transformation to peptide identifications

