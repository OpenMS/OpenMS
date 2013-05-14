from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DPosition cimport *
from DBoundingBox cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/ConvexHull2D.h>" namespace "OpenMS":

    cdef cppclass ConvexHull2D:

        ConvexHull2D() nogil except +
        ConvexHull2D(ConvexHull2D) nogil except + # wrap-ignore

        bool operator==(ConvexHull2D) nogil except +
        void clear()       nogil except +
        void compress()    nogil except +

        bool addPoint(DPosition2 point) nogil except +  # wrap-ignore
        void addPoints(libcpp_vector[DPosition2] points) nogil except +  # wrap-ignore

        void expandToBoundingBox() nogil except +
        bool encloses(DPosition2) nogil except + # wrap-ignore
        libcpp_vector[DPosition2] getHullPoints() nogil except + #wrap-ignore
        void setHullPoints(libcpp_vector[DPosition2] ) nogil except + #wrap-ignore
        DBoundingBox2 getBoundingBox() nogil except + #wrap-ignore

