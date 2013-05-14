
cdef extern from "<OpenMS/DATASTRUCTURES/CVMappings.h>" namespace "OpenMS":

    cdef cppclass CVMappings:

        CVMappings() nogil except +
        CVMappings(CVMappings) nogil except +

