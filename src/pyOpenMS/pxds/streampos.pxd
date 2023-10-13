
cdef extern from "<iostream>" namespace "std":
    cdef cppclass streampos:
        streampos() except + nogil 
        streampos(streampos &) except + nogil 
