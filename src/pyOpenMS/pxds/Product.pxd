from Types cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/METADATA/Product.h>" namespace "OpenMS":

    cdef cppclass Product:
        # wrap-doc:
        #  This class describes the product isolation window for special scan types, such as MRM

        Product() except + nogil 
        Product(Product &) except + nogil 

        bool operator==(Product) except + nogil 
        bool operator!=(Product) except + nogil 

        double getMZ() except + nogil  # wrap-doc:Returns the target m/z
        void setMZ(double) except + nogil  # wrap-doc:Sets the target m/z

        double getIsolationWindowLowerOffset() except + nogil  # wrap-doc:Returns the lower offset from the target m/z
        void  setIsolationWindowLowerOffset(double bound) except + nogil  # wrap-doc:Sets the lower offset from the target m/z

        double getIsolationWindowUpperOffset() except + nogil  # wrap-doc:Returns the upper offset from the target m/z
        void  setIsolationWindowUpperOffset(double bound) except + nogil  # wrap-doc:Sets the upper offset from the target m/z
