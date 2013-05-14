from InterpolationModel cimport *
from String cimport *
# from BasicStatistics cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>" namespace "OpenMS":
    
    cdef cppclass EmgModel(InterpolationModel):
        # wrap-inherits:
        #  InterpolationModel
        EmgModel() nogil except +
        EmgModel(EmgModel) nogil except +
        # BaseModel[ 1 ] * create() nogil except +
        String getProductName() nogil except +

        double getIntensity(DPosition1 &pos) nogil except +
        double getIntensity(double coord) nogil except +
        double getScalingFactor() nogil except +
        void setOffset(double offset) nogil except +
        double getCenter() nogil except +
        void setSamples() nogil except +
        void setInterpolationStep(double interpolation_step) nogil except +
        void setScalingFactor(double scaling) nogil except +
