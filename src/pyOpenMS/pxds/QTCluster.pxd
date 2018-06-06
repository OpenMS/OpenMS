from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from libcpp.set cimport set as libcpp_set
# from GridFeature cimport *
from AASequence cimport *

# typedef OpenMSBoost::unordered_map<Size, std::multimap<double, GridFeature*> > NeighborMap;

cdef extern from "<OpenMS/DATASTRUCTURES/QTCluster.h>" namespace "OpenMS":
    
    cdef cppclass QTCluster "OpenMS::QTCluster":
        QTCluster(QTCluster) nogil except + #wrap-ignore

        # POINTER #  QTCluster(GridFeature * center_point, Size num_maps, double max_distance, bool use_IDs) nogil except +
        double getCenterRT() nogil except +
        double getCenterMZ() nogil except +
        Int getXCoord() nogil except +
        Int getYCoord() nogil except +
        Size size() nogil except +
        bool operator<(QTCluster & cluster) nogil except +
        # POINTER # void add(GridFeature * element, double distance) nogil except +
        # POINTER # void getElements(libcpp_map[ Size, GridFeature * ] & elements) nogil except +
        # POINTER # bool update(libcpp_map[ Size, GridFeature * ] & removed) nogil except +
        double getQuality() nogil except +
        libcpp_set[ AASequence ]  getAnnotations() nogil except +
        void setInvalid() nogil except +
        bool isInvalid() nogil except +
        void initializeCluster() nogil except +
        void finalizeCluster() nogil except +
        # NAMESPACE # # POINTER # OpenMSBoost::unordered_map[ Size, libcpp_vector[ GridFeature * ] ] getAllNeighbors() nogil except +
        # NeighborMap getNeighbors() nogil except +

