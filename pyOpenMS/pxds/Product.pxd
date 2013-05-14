from Types cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/METADATA/Product.h>" namespace "OpenMS":

    cdef cppclass Product:
        

        Product()    nogil except +
        Product(Product)    nogil except +

        bool operator==(Product) nogil except +
        bool operator!=(Product) nogil except +

        DoubleReal getMZ()
        void setMZ(DoubleReal)

        DoubleReal getIsolationWindowLowerOffset()
        void  setIsolationWindowLowerOffset(DoubleReal bound)

        DoubleReal getIsolationWindowUpperOffset()
        void  setIsolationWindowUpperOffset(DoubleReal bound)
