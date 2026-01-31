from libcpp cimport *
from Types cimport *
from DPosition cimport *

cdef extern from "<OpenMS/KERNEL/Peak1D.h>" namespace "OpenMS":

    cdef cppclass Peak1D:
        # wrap-hash:
        #  std
        # wrap-doc:
        #  A 1-dimensional raw data point or peak.
        #  This data structure is intended for continuous data or peak data.
        #  If you want to annotate single peaks with meta data, use RichPeak1D instead.

        Peak1D() except + nogil  # wrap-doc:Default constructor
        Peak1D(Peak1D &) except + nogil  # wrap-doc:Copy constructor

        # We will not catch C++ exceptions for get/set methods for performance
        # reasons (no memory allocation is involved).
        float getIntensity() nogil  # wrap-doc:Returns the intensity (height) of the peak
        double getMZ() nogil  # wrap-doc:Returns the m/z (mass-to-charge) value of the peak
        void setMZ(double) nogil  # wrap-doc:Sets the m/z (mass-to-charge) value of the peak
        void setIntensity(float) nogil  # wrap-doc:Sets the intensity (height) of the peak
        bool operator==(Peak1D) except + nogil  # wrap-doc:Equality operator
        bool operator!=(Peak1D) except + nogil  # wrap-doc:Inequality operator
        double getPos() nogil  # wrap-doc:Returns the position (alias for getMZ)
        void setPos(double pos) nogil  # wrap-doc:Sets the position (alias for setMZ)
        # DPosition1 getPosition() except + nogil  # wrap-ignore
        # void setPosition(DPosition1 position) except + nogil  # wrap-ignore

