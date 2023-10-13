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
        QTCluster() except + nogil  # wrap-ignore
        QTCluster(QTCluster &) except + nogil 

        # POINTER #  QTCluster(GridFeature * center_point, Size num_maps, double max_distance, bool use_IDs) except + nogil 
        double getCenterRT() except + nogil  # wrap-doc:Returns the RT value of the cluster
        double getCenterMZ() except + nogil  # wrap-doc:Returns the m/z value of the cluster center
        Int getXCoord() except + nogil  # wrap-doc:Returns the x coordinate in the grid
        Int getYCoord() except + nogil  # wrap-doc:Returns the y coordinate in the grid
        Size size() except + nogil  # wrap-doc:Returns the size of the cluster (number of elements, incl. center)
        bool operator<(QTCluster & cluster) except + nogil 
        # POINTER # void add(GridFeature * element, double distance) except + nogil 
        # POINTER # void getElements(libcpp_map[ Size, GridFeature * ] & elements) except + nogil 
        # POINTER # bool update(libcpp_map[ Size, GridFeature * ] & removed) except + nogil 
        double getQuality() except + nogil  # wrap-doc:Returns the cluster quality and recomputes if necessary
        libcpp_set[ AASequence ]  getAnnotations() except + nogil  # wrap-doc:Returns the set of peptide sequences annotated to the cluster center
        void setInvalid() except + nogil  # wrap-doc:Sets current cluster as invalid (also frees some memory)
        bool isInvalid() except + nogil  # wrap-doc:Whether current cluster is invalid
        void initializeCluster() except + nogil  # wrap-doc:Has to be called before adding elements (calling QTCluster::add)
        void finalizeCluster() except + nogil  # wrap-doc:Has to be called after adding elements (after calling QTCluster::add one or multiple times)
        # NAMESPACE # # POINTER # OpenMSBoost::unordered_map[ Size, libcpp_vector[ GridFeature * ] ] getAllNeighbors() except + nogil 
        # NeighborMap getNeighbors() except + nogil 
