from Types cimport *
from libcpp cimport bool
# from DIntervalBase cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/DRange.h>" namespace "OpenMS":
    
    # TODO write fxn
    cdef cppclass DRange1 "OpenMS::DRange<1>":
        # wrap-ignore
        DRange1() nogil except +
        DRange1(DRange1) nogil except +
        ## DRange(PositionType & lower, PositionType & upper) nogil except +
        ## DRange(Base & range_) nogil except +
        ## DRange(CoordinateType minx, CoordinateType miny, CoordinateType maxx, CoordinateType maxy) nogil except +
        ## bool operator==(DRange & rhs) nogil except +
        ## bool operator==(Base & rhs) nogil except +
        ## bool encloses(PositionType & position) nogil except +
        ## bool encloses(CoordinateType x, CoordinateType y) nogil except +
        ## DRange united(DRange[ D ] & other_range) nogil except +
        ## DRangeIntersection intersects(DRange & range_) nogil except +
        ## bool isIntersected(DRange & range_) nogil except +
        ## bool isEmpty() nogil except +
        ## #  OpenMS::Internal::DIntervalBase< D >

cdef extern from "<OpenMS/DATASTRUCTURES/DRange.h>" namespace "OpenMS::DRange":
    cdef enum DRangeIntersection "OpenMS::DRange::DRangeIntersection":
        Disjoint
        Intersects
        Inside

