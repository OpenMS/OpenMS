from Types cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/METADATA/Product.h>" namespace "OpenMS":

    cdef cppclass Product:


        Product()    nogil except +
        Product(Product)    nogil except +

        bool operator==(Product) nogil except +
        bool operator!=(Product) nogil except +

        double getMZ() nogil except +
        void setMZ(double) nogil except +

        double getIsolationWindowLowerOffset() nogil except +
        void  setIsolationWindowLowerOffset(double bound) nogil except +

        double getIsolationWindowUpperOffset() nogil except +
        void  setIsolationWindowUpperOffset(double bound) nogil except +
