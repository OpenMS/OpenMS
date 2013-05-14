# from BaseModel cimport *
# from LinearInterpolation cimport *
from DPosition cimport *

ctypedef DoubleReal IntensityType
ctypedef DPosition1 PositionType
ctypedef DoubleReal CoordinateType
# ctypedef Math::LinearInterpolation<DoubleReal> LinearInterpolation

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>" namespace "OpenMS":
    
    cdef cppclass InterpolationModel "OpenMS::InterpolationModel":
        InterpolationModel() nogil except +
        InterpolationModel(InterpolationModel) nogil except +
        
        # InterpolationModel() nogil except +
        # InterpolationModel(InterpolationModel) nogil except +
        # double getIntensity(DPosition1 &pos) nogil except +
        # double getIntensity(double coord) nogil except +
        # double getScalingFactor() nogil except +
        # void setOffset(double offset) nogil except +
        # double getCenter() nogil except +
        # void setSamples() nogil except +
        # void setInterpolationStep(double interpolation_step) nogil except +
        # void setScalingFactor(double scaling) nogil except +

        # LinearInterpolation  getInterpolation() nogil except +
        # void getSamples(SamplesType &cont) nogil except +

