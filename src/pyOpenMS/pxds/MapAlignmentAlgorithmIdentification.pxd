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
        #   DefaultParamHandler
        #   ProgressLogger

        MapAlignmentAlgorithmIdentification() except + nogil 
        # private
        MapAlignmentAlgorithmIdentification(MapAlignmentAlgorithmIdentification &) except + nogil  # wrap-ignore
        
        void align(libcpp_vector[MSExperiment]&, libcpp_vector[TransformationDescription]&, int) except + nogil 
        void align(libcpp_vector[FeatureMap]&, libcpp_vector[TransformationDescription]&, int) except + nogil 
        void align(libcpp_vector[ConsensusMap]&, libcpp_vector[TransformationDescription]&, int) except + nogil 
        # TODO nested STL
        void align(libcpp_vector[libcpp_vector[PeptideIdentification]]& ids, libcpp_vector[TransformationDescription]& trafos, int ref_index) except + nogil  #wrap-ignore

        void setReference(MSExperiment&) except + nogil 
        void setReference(FeatureMap&) except + nogil 
        void setReference(ConsensusMap&) except + nogil 
        void setReference(libcpp_vector[PeptideIdentification]&) except + nogil 
