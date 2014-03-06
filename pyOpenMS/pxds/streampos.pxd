
cdef extern from "<iostream>" namespace "std":
    cdef cppclass streampos:
        streampos() nogil except +
        streampos(streampos &) nogil except +
