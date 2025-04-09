from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DPosition cimport *
from DBoundingBox cimport *

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/DATASTRUCTURES/ConvexHull2D.h>" namespace "OpenMS":

    cdef cppclass ConvexHull2D:

        ConvexHull2D() except + nogil 
        ConvexHull2D(ConvexHull2D &) except + nogil 

        bool operator==(ConvexHull2D) except + nogil 
        void clear()       except + nogil  # wrap-doc:Removes all points
        Size compress()    except + nogil  # wrap-doc:Allows to reduce the disk/memory footprint of a hull
        void expandToBoundingBox() except + nogil  # wrap-doc:Expand a convex hull to its bounding box.

        # as cython/autowrap have problems handling integer template args
        # we wrap following methods manually in ../addons/ConvexHull2D.pyx:
        # better/todo: implement converter in ../converters !

        bool addPoint(DPosition2 point) except + nogil  # wrap-doc:Adds a point to the hull if it is not already contained. Returns if the point was added. This will trigger recomputation of the outer hull points (thus points set with setHullPoints() will be lost)
        void addPoints(libcpp_vector[DPosition2] points) except + nogil  # wrap-doc:Adds points to the hull if it is not already contained. This will trigger recomputation of the outer hull points (thus points set with setHullPoints() will be lost)

        bool encloses(DPosition2) except + nogil  # wrap-doc:Returns if the `point` lies in the feature hull
        libcpp_vector[DPosition2] getHullPoints() except + nogil  # wrap-doc:Accessor for the outer points
        void setHullPoints(libcpp_vector[DPosition2] ) except + nogil  # wrap-doc:Accessor for the outer(!) points (no checking is performed if this is actually a convex hull) 
        DBoundingBox2 getBoundingBox() except + nogil  # wrap-doc:Returns the bounding box of the feature hull points
