from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from libcpp.set cimport set as libcpp_set
# from GridFeature cimport *
from AASequence cimport *

# typedef OpenMSBoost::unordered_map<Size, std::multimap<double, GridFeature*> > NeighborMap;

cdef extern from "<OpenMS/DATASTRUCTURES/QTCluster.h>" namespace "OpenMS":
    
    cdef cppclass QTCluster "OpenMS::QTCluster":
        # deleted
        QTCluster() nogil except + # wrap-ignore
        QTCluster(QTCluster &) nogil except +

        # POINTER #  QTCluster(GridFeature * center_point, Size num_maps, double max_distance, bool use_IDs) nogil except +
        double getCenterRT() nogil except + # wrap-doc:Returns the RT value of the cluster
        double getCenterMZ() nogil except + # wrap-doc:Returns the m/z value of the cluster center
        Int getXCoord() nogil except + # wrap-doc:Returns the x coordinate in the grid
        Int getYCoord() nogil except + # wrap-doc:Returns the y coordinate in the grid
        Size size() nogil except + # wrap-doc:Returns the size of the cluster (number of elements, incl. center)
        bool operator<(QTCluster & cluster) nogil except +
        # POINTER # void add(GridFeature * element, double distance) nogil except +
        # POINTER # void getElements(libcpp_map[ Size, GridFeature * ] & elements) nogil except +
        # POINTER # bool update(libcpp_map[ Size, GridFeature * ] & removed) nogil except +
        double getQuality() nogil except + # wrap-doc:Returns the cluster quality and recomputes if necessary
        libcpp_set[ AASequence ]  getAnnotations() nogil except + # wrap-doc:Returns the set of peptide sequences annotated to the cluster center
        void setInvalid() nogil except + # wrap-doc:Sets current cluster as invalid (also frees some memory)
        bool isInvalid() nogil except + # wrap-doc:Whether current cluster is invalid
        void initializeCluster() nogil except + # wrap-doc:Has to be called before adding elements (calling QTCluster::add)
        void finalizeCluster() nogil except + # wrap-doc:Has to be called after adding elements (after calling QTCluster::add one or multiple times)
        # NAMESPACE # # POINTER # OpenMSBoost::unordered_map[ Size, libcpp_vector[ GridFeature * ] ] getAllNeighbors() nogil except +
        # NeighborMap getNeighbors() nogil except +
