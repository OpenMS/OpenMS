from DPosition cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/DBoundingBox.h>" namespace "OpenMS":

    # this is the way to declare a class with a int template parameter:
    # as an alternative one could use a similar ctypedef, but these can
    # not be renamed when cipmorting, which might cause trouble !

    cdef cppclass DBoundingBox2 "OpenMS::DBoundingBox<2>":
        DBoundingBox2() nogil except +
        DBoundingBox2(DBoundingBox2) nogil except +
        DPosition2 minPosition() nogil except +
        DPosition2 maxPosition() nogil except +




