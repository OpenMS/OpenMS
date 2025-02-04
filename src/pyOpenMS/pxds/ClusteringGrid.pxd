from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DPosition cimport *

cdef extern from "<OpenMS/ML/CLUSTERING/ClusteringGrid.h>" namespace "OpenMS":
    
    cdef cppclass ClusteringGrid "OpenMS::ClusteringGrid":
        ClusteringGrid(libcpp_vector[ double ] & grid_spacing_x, libcpp_vector[ double ] & grid_spacing_y) except + nogil 
        ClusteringGrid(ClusteringGrid &) except + nogil  # compiler
        libcpp_vector[ double ] getGridSpacingX() except + nogil 
        libcpp_vector[ double ] getGridSpacingY() except + nogil 
        void addCluster(libcpp_pair[int,int] cell_index, int & cluster_index) except + nogil  # wrap-doc:Adds a cluster to this grid cell
        void removeCluster(libcpp_pair[int,int] cell_index, int & cluster_index) except + nogil  # wrap-doc:Removes a cluster from this grid cell and removes the cell if no other cluster left
        void removeAllClusters() except + nogil  # wrap-doc:Removes all clusters from this grid (and hence all cells)
        # NAMESPACE # std::list[ int ] getClusters(CellIndex & cell_index) except + nogil 
        libcpp_pair[int,int] getIndex(DPosition2 position) except + nogil 
        bool isNonEmptyCell(libcpp_pair[int,int] cell_index) except + nogil  # wrap-doc:Checks if there are clusters at this cell index
        int getCellCount() except + nogil  # wrap-doc:Returns number of grid cells occupied by one or more clusters
