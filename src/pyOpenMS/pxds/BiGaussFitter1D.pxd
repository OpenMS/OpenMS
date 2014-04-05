# from MaxLikeliFitter1D cimport *
# from BasicStatistics cimport *
from String cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussFitter1D.h>" namespace "OpenMS":
    
    cdef cppclass BiGaussFitter1D "OpenMS::BiGaussFitter1D":
        BiGaussFitter1D() nogil except +
        BiGaussFitter1D(BiGaussFitter1D) nogil except +
        # QualityType fit1d(RawDataArrayType &range, InterpolationModel *&model) nogil except +
        # Fitter1D * create() nogil except +
        String getProductName() nogil except +

