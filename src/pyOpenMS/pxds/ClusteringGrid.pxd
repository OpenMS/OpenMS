from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DPosition cimport *

cdef extern from "<OpenMS/COMPARISON/CLUSTERING/ClusteringGrid.h>" namespace "OpenMS":
    
    cdef cppclass ClusteringGrid "OpenMS::ClusteringGrid":
        ClusteringGrid(libcpp_vector[ double ] & grid_spacing_x, libcpp_vector[ double ] & grid_spacing_y) nogil except +
        ClusteringGrid(ClusteringGrid &) nogil except + # compiler
        libcpp_vector[ double ] getGridSpacingX() nogil except +
        libcpp_vector[ double ] getGridSpacingY() nogil except +
        void addCluster(libcpp_pair[int,int] cell_index, int & cluster_index) nogil except + # wrap-doc:Adds a cluster to this grid cell
        void removeCluster(libcpp_pair[int,int] cell_index, int & cluster_index) nogil except + # wrap-doc:Removes a cluster from this grid cell and removes the cell if no other cluster left
        void removeAllClusters() nogil except + # wrap-doc:Removes all clusters from this grid (and hence all cells)
        # NAMESPACE # std::list[ int ] getClusters(CellIndex & cell_index) nogil except +
        libcpp_pair[int,int] getIndex(DPosition2 position) nogil except +
        bool isNonEmptyCell(libcpp_pair[int,int] cell_index) nogil except + # wrap-doc:Checks if there are clusters at this cell index
        int getCellCount() nogil except + # wrap-doc:Returns number of grid cells occupied by one or more clusters
