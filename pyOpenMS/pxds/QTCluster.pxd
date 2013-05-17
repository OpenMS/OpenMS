from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from libcpp.set cimport set as libcpp_set
# from GridFeature cimport *
from AASequence cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/QTCluster.h>" namespace "OpenMS":
    
    cdef cppclass QTCluster "OpenMS::QTCluster":
        QTCluster(QTCluster) nogil except + #wrap-ignore
        # POINTER #  QTCluster(GridFeature * center_point, Size num_maps, DoubleReal max_distance, bool use_IDs) nogil except +
        DoubleReal getCenterRT() nogil except +
        DoubleReal getCenterMZ() nogil except +
        Size size() nogil except +
        bool operator<(QTCluster & cluster) nogil except +
        # POINTER # void add(GridFeature * element, DoubleReal distance) nogil except +
        # POINTER # void getElements(libcpp_map[ Size, GridFeature * ] & elements) nogil except +
        # POINTER # bool update(libcpp_map[ Size, GridFeature * ] & removed) nogil except +
        DoubleReal getQuality() nogil except +
        libcpp_set[ AASequence ]  getAnnotations() nogil except +

