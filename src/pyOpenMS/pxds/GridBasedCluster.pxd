from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DPosition cimport *
from DBoundingBox cimport *


# typedef DPosition<2> Point;
# typedef DBoundingBox<2> Rectangle;

cdef extern from "<OpenMS/COMPARISON/CLUSTERING/GridBasedCluster.h>" namespace "OpenMS":
    
    cdef cppclass GridBasedCluster "OpenMS::GridBasedCluster":

        GridBasedCluster(DPosition2 centre,
            DBoundingBox2 bounding_box, 
            libcpp_vector[ int ] point_indices,
            int property_A, 

            libcpp_vector[ int ] properties_B) nogil except +

        GridBasedCluster(DPosition2 centre, 
           DBoundingBox2 bounding_box,
           libcpp_vector[ int ] point_indices) nogil except +

        GridBasedCluster(GridBasedCluster &) nogil except + # compiler

        DPosition2 getCentre() nogil except + # wrap-doc:Returns cluster centre
        DBoundingBox2 getBoundingBox() nogil except + # wrap-doc:Returns bounding box
        libcpp_vector[ int ] getPoints() nogil except + # wrap-doc:Returns indices of points in cluster
        int getPropertyA() nogil except + # wrap-doc:Returns property A
        libcpp_vector[ int ] getPropertiesB() nogil except + # wrap-doc:Returns properties B of all points


        # bool operator<(GridBasedCluster other) nogil except +
        # bool operator>(GridBasedCluster other) nogil except +
        # bool operator==(GridBasedCluster other) nogil except +
