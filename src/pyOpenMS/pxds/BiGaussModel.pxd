# from InterpolationModel cimport *
# from BasicStatistics cimport *
from String cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussModel.h>" namespace "OpenMS":
    
    cdef cppclass BiGaussModel "OpenMS::BiGaussModel":
        BiGaussModel() nogil except +
        BiGaussModel(BiGaussModel) nogil except + # wrap-ignore
        void setOffset(double offset) nogil except +
        void setSamples() nogil except +
        double getCenter() nogil except +
        # BaseModel[ 1 ] * create() nogil except +
        String getProductName() nogil except +

