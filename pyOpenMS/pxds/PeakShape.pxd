from Types cimport *
from libcpp cimport bool
from Types cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>" namespace "OpenMS":
    
    cdef cppclass PeakShape "OpenMS::PeakShape":
        PeakShape() nogil except +
        PeakShape(PeakShape) nogil except +
        DoubleReal height
        DoubleReal mz_position
        DoubleReal left_width
        DoubleReal right_width
        DoubleReal area
        DoubleReal r_value
        DoubleReal signal_to_noise
        PeakShape_Type type
        # PeakShape(DoubleReal height_, DoubleReal mz_position_, DoubleReal left_width_, DoubleReal right_width_, DoubleReal area_, PeakIterator left_, PeakIterator right_, PeakShape_Type type_) nogil except +
        PeakShape(DoubleReal height_, DoubleReal mz_position_, DoubleReal left_width_, DoubleReal right_width_, DoubleReal area_, PeakShape_Type type_) nogil except +
        bool operator==(PeakShape & rhs) nogil except +
        bool operator!=(PeakShape & rhs) nogil except +
        # DoubleReal operator()(DoubleReal x) nogil except +
        DoubleReal getSymmetricMeasure() nogil except +
        DoubleReal getFWHM() nogil except +
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

