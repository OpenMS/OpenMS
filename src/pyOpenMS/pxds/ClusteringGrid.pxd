from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DPosition cimport *

cdef extern from "<OpenMS/COMPARISON/CLUSTERING/ClusteringGrid.h>" namespace "OpenMS":
    
    cdef cppclass ClusteringGrid "OpenMS::ClusteringGrid":
        ClusteringGrid(ClusteringGrid) nogil except + #wrap-ignore
        ClusteringGrid(libcpp_vector[ double ] & grid_spacing_x, libcpp_vector[ double ] & grid_spacing_y) nogil except +
        libcpp_vector[ double ] getGridSpacingX() nogil except +
        libcpp_vector[ double ] getGridSpacingY() nogil except +
        void addCluster(libcpp_pair[int,int] cell_index, int & cluster_index) nogil except +
        void removeCluster(libcpp_pair[int,int] cell_index, int & cluster_index) nogil except +
        void removeAllClusters() nogil except +
        # NAMESPACE # std::list[ int ] getClusters(CellIndex & cell_index) nogil except +
        libcpp_pair[int,int] getIndex(DPosition2 position) nogil except +
        bool isNonEmptyCell(libcpp_pair[int,int] cell_index) nogil except +
        int getCellCount() nogil except +
