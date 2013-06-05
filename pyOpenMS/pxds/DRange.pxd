from Types cimport *
from DPosition cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/DATASTRUCTURES/DRange.h>" namespace "OpenMS":
    
    cdef cppclass DRange1 "OpenMS::DRange<1>":
        DRange1() nogil except +
        DRange1(DRange1) nogil except +
        DRange1(DPosition1 lower, DPosition1 upper) nogil except +
        bool operator==(DRange1 & rhs) nogil except +
        bool encloses(DPosition1 & position) nogil except +
        DRange1 united(DRange1 other_range) nogil except +
        # DRangeIntersection intersects(DRange1 & range_) nogil except +
        bool isIntersected(DRange1 & range_) nogil except +
        bool isEmpty() nogil except +

    cdef cppclass DRange2 "OpenMS::DRange<2>":
        DRange2() nogil except +
        DRange2(DRange2) nogil except +
        DRange2(DPosition2 lower, DPosition2 upper) nogil except +
        DRange2(DoubleReal minx, DoubleReal miny, DoubleReal maxx, DoubleReal maxy) nogil except +
        bool operator==(DRange2 & rhs) nogil except +
        # bool encloses(DPosition2 & position) nogil except +
        # bool encloses(DoubleReal x, DoubleReal y) nogil except +
        DRange2 united(DRange2 other_range) nogil except +
        # DRangeIntersection intersects(DRange2 & range_) nogil except +
        bool isIntersected(DRange2 & range_) nogil except +
        bool isEmpty() nogil except +

cdef extern from "<OpenMS/DATASTRUCTURES/DRange.h>" namespace "OpenMS::DRange":
    cdef enum DRangeIntersection "OpenMS::DRange::DRangeIntersection":
        Disjoint
        Intersects
        Inside

