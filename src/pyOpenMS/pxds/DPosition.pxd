from libcpp cimport bool
from Types cimport *

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/DATASTRUCTURES/DPosition.h>" namespace "OpenMS":

    # this is the way to declare a class with a int template parameter.

    # it is a bit laborous as one has to declare all methods twice,
    # but deriving from a common base class and giving an extra string
    # with explicit c++ type declarion does not work !

    # as an alternative one could use a similar ctypedef, but these can
    # not be renamed when cipmorting, which might cause trouble !

    cdef cppclass DPosition1 "OpenMS::DPosition<1> ":
        # wrap-ignore
        DPosition1()  nogil except +
        DPosition1(double)  nogil except +
        DPosition1(DPosition1)  nogil except +

        double & operator[](Size index) nogil except +
        bool operator==(DPosition1) nogil except +
        bool operator!=(DPosition1) nogil except +

        bool operator<(DPosition1) nogil except +
        bool operator<=(DPosition1) nogil except +

        bool operator>(DPosition1) nogil except +
        bool operator>=(DPosition1) nogil except +

        void clear() nogil except +
        Size size() nogil except +

    cdef cppclass DPosition2 "OpenMS::DPosition<2> ":
        # wrap-ignore
        DPosition2()  nogil except +
        DPosition2(DPosition2)  nogil except +
        DPosition2(double)  nogil except +
        DPosition2(double, double)  nogil except +

        double & operator[](Size index) nogil except +
        bool operator==(DPosition2) nogil except +
        bool operator!=(DPosition2) nogil except +

        bool operator<(DPosition2) nogil except +
        bool operator<=(DPosition2) nogil except +

        bool operator>(DPosition2) nogil except +
        bool operator>=(DPosition2) nogil except +

        void clear() nogil except +
        Size size() nogil except +
