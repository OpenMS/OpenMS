from Types cimport *
from KDTreeFeatureMaps cimport *
from TransformationModelLowess cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmKD.h>" namespace "OpenMS":
    
    cdef cppclass MapAlignmentAlgorithmKD "OpenMS::MapAlignmentAlgorithmKD":
        # private
        MapAlignmentAlgorithmKD() nogil except + # wrap-ignore

        MapAlignmentAlgorithmKD(MapAlignmentAlgorithmKD) nogil except + # compiler

        MapAlignmentAlgorithmKD(Size num_maps, Param & param) nogil except +

        void addRTFitData(KDTreeFeatureMaps & kd_data) nogil except +

        void fitLOWESS() nogil except +

        void transform(KDTreeFeatureMaps & kd_data) nogil except +

