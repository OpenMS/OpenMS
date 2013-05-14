
# enum support in cython is tricky, as cython has problems with
# nested structrues / classes / ....
# so: set 'namespace' and 'ns' as seen below !

cdef extern from "<OpenMS/METADATA/IonSource.h>" namespace "OpenMS::IonSource":

     cdef enum Polarity:
            POLNULL, POSITIVE, NEGATIVE, SIZE_OF_POLARITY
