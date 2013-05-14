

cdef extern from "<OpenMS/METADATA/IonSource.h>" namespace "OpenMS::IonSource":

     cdef enum Polarity:
            POLNULL, POSITIVE, NEGATIVE, SIZE_OF_POLARITY
