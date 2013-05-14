from libcpp cimport bool
from Types cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/DPosition.h>" namespace "OpenMS":

    cdef cppclass DPosition[N]:
        # wrap-ignore

        DPosition()  nogil except +
        DPosition(DPosition[N])  nogil except +

        DoubleReal operator[](Size index)
        bool operator==(DPosition[N])
        bool operator!=(DPosition[N])

        bool operator<(DPosition[N])
        bool operator<=(DPosition[N])

        bool operator>(DPosition[N])
        bool operator>=(DPosition[N])

        void clear()
        Size size()

    #ctypedef DPosition[unsigned] DPosition1 "DPosition<1>"
