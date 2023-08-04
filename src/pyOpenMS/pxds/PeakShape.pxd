from Types cimport *
from libcpp cimport bool
from Types cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>" namespace "OpenMS":
    
    cdef cppclass PeakShape "OpenMS::PeakShape":
        # wrap-doc:
            #  Internal representation of a peak shape (used by the PeakPickerCWT)
            #  
            #  It defines an asymmetric Lorentzian and asymmetric hyperbolic squared secan function

        PeakShape() except + nogil 
        PeakShape(PeakShape &) except + nogil 
        double height
        double mz_position
        double left_width
        double right_width
        double area
        double r_value
        double signal_to_noise
        PeakShape_Type type
        # PeakShape(double height_, double mz_position_, double left_width_, double right_width_, double area_, PeakIterator left_, PeakIterator right_, PeakShape_Type type_) except + nogil 
        PeakShape(double height_, double mz_position_, double left_width_, double right_width_, double area_, PeakShape_Type type_) except + nogil 
        bool operator==(PeakShape & rhs) except + nogil 
        bool operator!=(PeakShape & rhs) except + nogil 
        # double operator()(double x) except + nogil 
        double getSymmetricMeasure() except + nogil  # wrap-doc:Computes symmetry measure of the peak shape, which is corresponds to the ratio of the left and right width parameters
        double getFWHM() except + nogil  # wrap-doc:Estimates the full width at half maximum
        bool iteratorsSet() except + nogil  # wrap-doc:Check if endpoint iterators are provided
        # PeakIterator getLeftEndpoint() except + nogil 
        # void setLeftEndpoint(PeakIterator left_endpoint) except + nogil 
        # PeakIterator getRightEndpoint() except + nogil 
        # void setRightEndpoint(PeakIterator right_endpoint) except + nogil 

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>" namespace "OpenMS::PeakShape":
    cdef enum PeakShape_Type "OpenMS::PeakShape::Type":
        #wrap-attach:
        #   PeakShape
        LORENTZ_PEAK
        SECH_PEAK
        UNDEFINED

