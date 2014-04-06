from Types cimport *
from libcpp cimport bool
from Types cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>" namespace "OpenMS":
    
    cdef cppclass PeakShape "OpenMS::PeakShape":
        PeakShape() nogil except +
        PeakShape(PeakShape) nogil except +
        double height
        double mz_position
        double left_width
        double right_width
        double area
        double r_value
        double signal_to_noise
        PeakShape_Type type
        # PeakShape(double height_, double mz_position_, double left_width_, double right_width_, double area_, PeakIterator left_, PeakIterator right_, PeakShape_Type type_) nogil except +
        PeakShape(double height_, double mz_position_, double left_width_, double right_width_, double area_, PeakShape_Type type_) nogil except +
        bool operator==(PeakShape & rhs) nogil except +
        bool operator!=(PeakShape & rhs) nogil except +
        # double operator()(double x) nogil except +
        double getSymmetricMeasure() nogil except +
        double getFWHM() nogil except +
        bool iteratorsSet() nogil except +
        # PeakIterator getLeftEndpoint() nogil except +
        # void setLeftEndpoint(PeakIterator left_endpoint) nogil except +
        # PeakIterator getRightEndpoint() nogil except +
        # void setRightEndpoint(PeakIterator right_endpoint) nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>" namespace "OpenMS::PeakShape":
    cdef enum PeakShape_Type "OpenMS::PeakShape::Type":
        #wrap-attach:
        #    PeakShape
        LORENTZ_PEAK
        SECH_PEAK
        UNDEFINED

