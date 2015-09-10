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

        void transformRetentionTimes(MSExperiment[Peak1D,ChromatogramPeak]&, TransformationDescription&, bool) nogil except +

        void transformRetentionTimes(FeatureMap&, TransformationDescription&, bool) nogil except +

        void transformRetentionTimes(ConsensusMap&, TransformationDescription&, bool) nogil except +

        void transformRetentionTimes(libcpp_vector[PeptideIdentification]&, TransformationDescription&, bool) nogil except +

