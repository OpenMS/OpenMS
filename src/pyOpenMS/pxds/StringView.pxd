from Types cimport *
from String cimport *
from libcpp.string cimport string as libcpp_string
from libcpp cimport bool

cdef extern from "<OpenMS/DATASTRUCTURES/String.h>" namespace "OpenMS":
    cdef cppclass StringView:

        StringView() except + nogil  # TODO
        StringView(const libcpp_string &) except + nogil 
        StringView(StringView &) except + nogil 

        bool operator<(StringView other) except + nogil 
        StringView substr(Size start, Size end)  except + nogil 
        Size size()  except + nogil 
        String getString()  except + nogil 
