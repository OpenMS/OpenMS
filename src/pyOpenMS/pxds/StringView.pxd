from Types cimport *
from String cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/DATASTRUCTURES/String.h>" namespace "OpenMS":
    cdef cppclass StringView:

        StringView() nogil except +
        StringView(const String& s) nogil except +
        StringView(const StringView& s) nogil except +

        bool operator<(const StringView other) nogil except +
        StringView substr(Size start, Size end)  nogil except +
        Size size()  nogil except +
        String getString()  nogil except +
