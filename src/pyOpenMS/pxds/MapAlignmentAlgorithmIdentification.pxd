from libcpp.vector cimport vector as libcpp_vector

from ChromatogramPeak cimport *
from ConsensusMap cimport *
from DefaultParamHandler cimport *
from Feature cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from Param cimport *
from Peak1D cimport *
from PeptideIdentification cimport *
from ProgressLogger cimport *
from TransformationDescription cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>" namespace "OpenMS":

    cdef cppclass MapAlignmentAlgorithmIdentification(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        MapAlignmentAlgorithmIdentification() nogil except +
        
        void align(libcpp_vector[MSExperiment[Peak1D,ChromatogramPeak]]&, libcpp_vector[TransformationDescription]&, int) nogil except +
        void align(libcpp_vector[FeatureMap]&, libcpp_vector[TransformationDescription]&, int) nogil except +
        void align(libcpp_vector[ConsensusMap]&, libcpp_vector[TransformationDescription]&, int) nogil except +
        # TODO nested STL
        void align(libcpp_vector[libcpp_vector[PeptideIdentification]]& ids, libcpp_vector[TransformationDescription]& trafos, int ref_index) nogil except + #wrap-ignore

        void setReference(MSExperiment[Peak1D,ChromatogramPeak]&) nogil except +
        void setReference(FeatureMap&) nogil except +
        void setReference(ConsensusMap&) nogil except +
        void setReference(libcpp_vector[PeptideIdentification]&) nogil except +
