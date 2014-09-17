from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DPosition cimport *
from DBoundingBox cimport *


# typedef DPosition<2> Point;
# typedef DBoundingBox<2> Rectangle;

cdef extern from "<OpenMS/COMPARISON/CLUSTERING/GridBasedCluster.h>" namespace "OpenMS":
    
    cdef cppclass GridBasedCluster "OpenMS::GridBasedCluster":

        GridBasedCluster(GridBasedCluster) nogil except + #wrap-ignore
        GridBasedCluster(DPosition2 centre,
            DBoundingBox2 bounding_box, 
            libcpp_vector[ int ] point_indices,
            int property_A, 
            libcpp_vector[ int ] properties_B) nogil except +
        GridBasedCluster(DPosition2 centre, 
           DBoundingBox2 bounding_box,
           libcpp_vector[ int ] point_indices) nogil except +

        DPosition2 getCentre() nogil except +
        DBoundingBox2 getBoundingBox() nogil except +
        libcpp_vector[ int ] getPoints() nogil except +
        int getPropertyA() nogil except +
        libcpp_vector[ int ] getPropertiesB() nogil except +

        # bool operator<(GridBasedCluster other) nogil except +
        # bool operator>(GridBasedCluster other) nogil except +
        # bool operator==(GridBasedCluster other) nogil except +


