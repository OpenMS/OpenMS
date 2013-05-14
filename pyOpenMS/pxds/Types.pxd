from libc.stddef cimport *
from ctime cimport *

cdef extern from "<OpenMS/CONCEPT/Types.h>" namespace "OpenMS":

    ctypedef signed int Int32
    ctypedef signed long Int64
    ctypedef unsigned int UInt32
    ctypedef unsigned long UInt64
    ctypedef time_t   Time
    ctypedef unsigned int UInt
    ctypedef int      Int
    ctypedef float    Real
    ctypedef double   DoubleReal
    ctypedef unsigned long   UID
    ctypedef size_t    Size
    ctypedef ptrdiff_t SignedSize
