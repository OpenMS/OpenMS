from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DPosition cimport *
from DBoundingBox cimport *


# typedef DPosition<2> Point;
# typedef DBoundingBox<2> Rectangle;

cdef extern from "<OpenMS/ML/CLUSTERING/GridBasedCluster.h>" namespace "OpenMS":
    
    cdef cppclass GridBasedCluster "OpenMS::GridBasedCluster":

        GridBasedCluster(DPosition2 centre,
            DBoundingBox2 bounding_box, 
            libcpp_vector[ int ] point_indices,
            int property_A, 

            libcpp_vector[ int ] properties_B) except + nogil 

        GridBasedCluster(DPosition2 centre, 
           DBoundingBox2 bounding_box,
           libcpp_vector[ int ] point_indices) except + nogil 

        GridBasedCluster(GridBasedCluster &) except + nogil  # compiler

        DPosition2 getCentre() except + nogil  # wrap-doc:Returns cluster centre
        DBoundingBox2 getBoundingBox() except + nogil  # wrap-doc:Returns bounding box
        libcpp_vector[ int ] getPoints() except + nogil  # wrap-doc:Returns indices of points in cluster
        int getPropertyA() except + nogil  # wrap-doc:Returns property A
        libcpp_vector[ int ] getPropertiesB() except + nogil  # wrap-doc:Returns properties B of all points


        # bool operator<(GridBasedCluster other) except + nogil 
        # bool operator>(GridBasedCluster other) except + nogil 
        # bool operator==(GridBasedCluster other) except + nogil 
