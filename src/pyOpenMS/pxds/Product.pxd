from Types cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/METADATA/Product.h>" namespace "OpenMS":

    cdef cppclass Product:
        # wrap-doc:
        #   This class describes the product isolation window for special scan types, such as MRM

        Product() nogil except +
        Product(Product &) nogil except +

        bool operator==(Product) nogil except +
        bool operator!=(Product) nogil except +

        double getMZ() nogil except + # wrap-doc:Returns the target m/z
        void setMZ(double) nogil except + # wrap-doc:Sets the target m/z

        double getIsolationWindowLowerOffset() nogil except + # wrap-doc:Returns the lower offset from the target m/z
        void  setIsolationWindowLowerOffset(double bound) nogil except + # wrap-doc:Sets the lower offset from the target m/z

        double getIsolationWindowUpperOffset() nogil except + # wrap-doc:Returns the upper offset from the target m/z
        void  setIsolationWindowUpperOffset(double bound) nogil except + # wrap-doc:Sets the upper offset from the target m/z
