# from MaxLikeliFitter1D cimport *
# from BasicStatistics cimport *
from String cimport *

cdef extern from "<OpenMS/FEATUREFINDER/BiGaussFitter1D.h>" namespace "OpenMS":
    
    cdef cppclass BiGaussFitter1D "OpenMS::BiGaussFitter1D":
        # wrap-doc:
        #  BiGaussian distribution fitter (1-dim.) approximated using linear
        #  interpolation
        BiGaussFitter1D() except + nogil 
        BiGaussFitter1D(BiGaussFitter1D &) except + nogil 
        # QualityType fit1d(RawDataArrayType &range, InterpolationModel *&model) except + nogil 
        # Fitter1D * create() except + nogil 
       
