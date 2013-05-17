from Types cimport *
from DPosition cimport *
from libcpp cimport bool
# from DIntervalBase cimport *

# CoordinateType DoubleReal
# PositionType DPosition<D>
cdef extern from "<OpenMS/DATASTRUCTURES/DRange.h>" namespace "OpenMS":
    
    # TODO write fxn
    cdef cppclass DRange1 "OpenMS::DRange<1>":
        DRange1() nogil except +
        DRange1(DRange1) nogil except +
        # DRange(DPosition1 & lower, DPosition1 & upper) nogil except +
        ## DRange(Base & range_) nogil except +
        ## DRange(CoordinateType minx, CoordinateType miny, CoordinateType maxx, CoordinateType maxy) nogil except +
        bool operator==(DRange1 & rhs) nogil except +
        ## bool operator==(Base & rhs) nogil except +
        bool encloses(DPosition1 & position) nogil except +
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

