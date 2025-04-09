from Types cimport *
from DPosition cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/DATASTRUCTURES/DRange.h>" namespace "OpenMS":
    
    cdef cppclass DRange1 "OpenMS::DRange<1> ":
        DRange1() except + nogil  # TODO
        DRange1(DRange1 &) except + nogil 
        DRange1(DPosition1 lower, DPosition1 upper) except + nogil 
        bool operator==(DRange1 & rhs) except + nogil 
        bool encloses(DPosition1 & position) except + nogil 
        DRange1 united(DRange1 other_range) except + nogil 
        # DRangeIntersection intersects(DRange1 & range_) except + nogil 
        bool isIntersected(DRange1 & range_) except + nogil 
        bool isEmpty() except + nogil 

    cdef cppclass DRange2 "OpenMS::DRange<2> ":
        DRange2() except + nogil  # TODO
        DRange2(DRange2 &) except + nogil 
        DRange2(DPosition2 lower, DPosition2 upper) except + nogil 
        DRange2(double minx, double miny, double maxx, double maxy) except + nogil 
        bool operator==(DRange2 & rhs) except + nogil 
        # bool encloses(DPosition2 & position) except + nogil 
        # bool encloses(double x, double y) except + nogil 
        DRange2 united(DRange2 other_range) except + nogil 
        # DRangeIntersection intersects(DRange2 & range_) except + nogil 
        bool isIntersected(DRange2 & range_) except + nogil 
        bool isEmpty() except + nogil 

cdef extern from "<OpenMS/DATASTRUCTURES/DRange.h>" namespace "OpenMS::DRange":
    cdef enum DRangeIntersection "OpenMS::DRange::DRangeIntersection":
        Disjoint
        Intersects
        Inside

