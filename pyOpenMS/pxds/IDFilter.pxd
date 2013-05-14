from String cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/FILTERING/ID/IDFilter.h>" namespace "OpenMS":

    cdef cppclass IDFilter:

        IDFilter()                  nogil except +
        IDFilter(IDFilter)   nogil except + #wrap-ignore

        # TODO 

