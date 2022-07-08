from Types cimport *
from KDTreeFeatureMaps cimport *
from TransformationModelLowess cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmKD.h>" namespace "OpenMS":
    
    cdef cppclass MapAlignmentAlgorithmKD "OpenMS::MapAlignmentAlgorithmKD":
        # wrap-doc:
            #   An efficient reference-free feature map alignment algorithm for unlabeled data
            #   -----
            #   This algorithm uses a kd-tree to efficiently compute conflict-free connected components (CCC)
            #   in a compatibility graph on feature data. This graph is comprised of nodes corresponding
            #   to features and edges connecting features f and f' iff both are within each other's tolerance
            #   windows (wrt. RT and m/z difference). CCCs are those CCs that do not contain multiple features
            #   from the same input map, and whose features all have the same charge state
            #   -----
            #   All CCCs above a user-specified minimum size are considered true sets of corresponding features
            #   and based on these, LOWESS transformations are computed for each input map such that the average
            #   deviation from the mean retention time within all CCCs is minimized

        # private
        MapAlignmentAlgorithmKD() nogil except + # wrap-ignore
        MapAlignmentAlgorithmKD(Size num_maps, Param & param) nogil except +

        MapAlignmentAlgorithmKD(MapAlignmentAlgorithmKD &) nogil except + # compiler

        void addRTFitData(KDTreeFeatureMaps & kd_data) nogil except + # wrap-doc:Compute data points needed for RT transformation in the current `kd_data`, add to `fit_data_`
        void fitLOWESS() nogil except + # wrap-doc:Fit LOWESS to fit_data_, store final models in `transformations_`
        void transform(KDTreeFeatureMaps & kd_data) nogil except + # wrap-doc:Transform RTs for `kd_data`
