from String cimport *

cdef extern from "<OpenMS/FORMAT/TXTToMzMLConverter.h>" namespace "OpenMS":

    cdef cppclass TXTToMzMLConverter:

        TXTToMzMLConverter() nogil except +
        TXTToMzMLConverter(TXTToMzMLConverter) nogil except +
