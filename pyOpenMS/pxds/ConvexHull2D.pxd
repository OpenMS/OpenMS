from libcpp cimport bool

cdef extern from "<OpenMS/DATASTRUCTURES/ConvexHull2D.h>" namespace "OpenMS":

    cdef cppclass ConvexHull2D:

        ConvexHull2D()
        ConvexHull2D(ConvexHull2D) # wrap-ignore

        bool operator==(ConvexHull2D)
        void clear()

