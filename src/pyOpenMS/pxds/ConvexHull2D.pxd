from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DPosition cimport *
from DBoundingBox cimport *

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/DATASTRUCTURES/ConvexHull2D.h>" namespace "OpenMS":

    cdef cppclass ConvexHull2D:

        ConvexHull2D() nogil except +
        ConvexHull2D(ConvexHull2D &) nogil except +

        bool operator==(ConvexHull2D) nogil except +
        void clear()       nogil except + # wrap-doc:Removes all points
        Size compress()    nogil except + # wrap-doc:Allows to reduce the disk/memory footprint of a hull
        void expandToBoundingBox() nogil except + # wrap-doc:Expand a convex hull to its bounding box.

        # as cython/autowrap have problems handling integer template args
        # we wrap following methods manually in ../addons/ConvexHull2D.pyx:
        # better/todo: implement converter in ../converters !

        bool addPoint(DPosition2 point) nogil except + # wrap-doc:Adds a point to the hull if it is not already contained. Returns if the point was added. This will trigger recomputation of the outer hull points (thus points set with setHullPoints() will be lost)
        void addPoints(libcpp_vector[DPosition2] points) nogil except + # wrap-doc:Adds points to the hull if it is not already contained. This will trigger recomputation of the outer hull points (thus points set with setHullPoints() will be lost)

        bool encloses(DPosition2) nogil except + # wrap-doc:Returns if the `point` lies in the feature hull
        libcpp_vector[DPosition2] getHullPoints() nogil except + # wrap-doc:Accessor for the outer points
        void setHullPoints(libcpp_vector[DPosition2] ) nogil except + # wrap-doc:Accessor for the outer(!) points (no checking is performed if this is actually a convex hull) 
        DBoundingBox2 getBoundingBox() nogil except + # wrap-doc:Returns the bounding box of the feature hull points
