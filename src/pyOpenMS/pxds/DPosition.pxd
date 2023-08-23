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
        DPosition1()  except + nogil  # TODO
        DPosition1(double)  except + nogil 
        DPosition1(DPosition1 &)  except + nogil 

        double & operator[](Size index) except + nogil 
        bool operator==(DPosition1) except + nogil 
        bool operator!=(DPosition1) except + nogil 

        bool operator<(DPosition1) except + nogil 
        bool operator<=(DPosition1) except + nogil 

        bool operator>(DPosition1) except + nogil 
        bool operator>=(DPosition1) except + nogil 

        void clear() except + nogil 
        Size size() except + nogil 

    cdef cppclass DPosition2 "OpenMS::DPosition<2> ":
        # wrap-ignore
        DPosition2()  except + nogil  # TODO
        DPosition2(DPosition2 &)  except + nogil 
        DPosition2(double)  except + nogil 
        DPosition2(double, double)  except + nogil 

        double & operator[](Size index) except + nogil 
        bool operator==(DPosition2) except + nogil 
        bool operator!=(DPosition2) except + nogil 

        bool operator<(DPosition2) except + nogil 
        bool operator<=(DPosition2) except + nogil 

        bool operator>(DPosition2) except + nogil 
        bool operator>=(DPosition2) except + nogil 

        void clear() except + nogil 
        Size size() except + nogil 
