from libcpp cimport *
from Types cimport *
from DPosition cimport *

cdef extern from "<OpenMS/KERNEL/Peak2D.h>" namespace "OpenMS":

    cdef cppclass Peak2D:
        # wrap-doc:
            #   A 2-dimensional raw data point or peak.
            #   -----
            #   This data structure is intended for continuous data or peak data.
            #   If you want to annotated single peaks with meta data, use RichPeak2D instead

        Peak2D() nogil except +
        Peak2D(Peak2D &) nogil except +
        #Peak2D(DPosition2 &, float) nogil except +

        float getIntensity() nogil except + # wrap-doc:Returns the data point intensity (height)
        double getMZ() nogil except + # wrap-doc:Returns the m/z coordinate (index 1)
        double getRT() nogil except + # wrap-doc:Returns the RT coordinate (index 0)
        void setMZ(double) nogil except + # wrap-doc:Returns the m/z coordinate (index 1)
        void setRT(double) nogil except + # wrap-doc:Returns the RT coordinate (index 0)
        void setIntensity(float) nogil except + # wrap-doc:Returns the data point intensity (height)
        bool operator==(Peak2D) nogil except +
        bool operator!=(Peak2D) nogil except +
        # DPosition2 getPosition()
        # void setPosition(DPosition2 & position)

cdef extern from "<OpenMS/KERNEL/Peak2D.h>" namespace "OpenMS::Peak2D":
    
    cdef enum DimensionDescription "OpenMS::Peak2D::DimensionDescription":
        RT
        MZ
        DIMENSION
