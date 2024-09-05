# from InterpolationModel cimport *
# from BasicStatistics cimport *
from String cimport *

cdef extern from "<OpenMS/FEATUREFINDER/BiGaussModel.h>" namespace "OpenMS":
    
    cdef cppclass BiGaussModel "OpenMS::BiGaussModel":
        BiGaussModel() except + nogil 
        BiGaussModel(BiGaussModel &) except + nogil 
        void setOffset(double offset) except + nogil 
        void setSamples() except + nogil 
        double getCenter() except + nogil 
        # BaseModel * create() except + nogil 
       