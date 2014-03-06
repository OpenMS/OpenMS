from Types cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/METADATA/Product.h>" namespace "OpenMS":

    cdef cppclass Product:


        Product()    nogil except +
        Product(Product)    nogil except +

        bool operator==(Product) nogil except +
        bool operator!=(Product) nogil except +

        DoubleReal getMZ() nogil except +
        void setMZ(DoubleReal) nogil except +

        DoubleReal getIsolationWindowLowerOffset() nogil except +
        void  setIsolationWindowLowerOffset(DoubleReal bound) nogil except +

        DoubleReal getIsolationWindowUpperOffset() nogil except +
        void  setIsolationWindowUpperOffset(DoubleReal bound) nogil except +
