from libc.stddef cimport *
from libc.stdint cimport *
from ctime cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
from libcpp.map cimport map as libcpp_map
from libcpp.set cimport set as libcpp_set
from libcpp.string cimport string as libcpp_string
from smart_ptr cimport shared_ptr

# include macros
cdef extern from "<OpenMS/config.h>":
    pass

cdef extern from "<OpenMS/CONCEPT/Types.h>" namespace "OpenMS":

    ctypedef int32_t Int32     "OPENMS_INT32_TYPE"
    ctypedef int64_t Int64    "OPENMS_INT64_TYPE"
    ctypedef uint32_t UInt32  "OPENMS_UINT32_TYPE"
    ctypedef uint64_t UInt64 "OPENMS_UINT64_TYPE"
    ctypedef time_t   Time
    ctypedef unsigned int UInt
    ctypedef int      Int
    ctypedef uint64_t UID "OPENMS_UINT64_TYPE"
    ctypedef size_t    Size
    ctypedef ptrdiff_t SignedSize
