# from MaxLikeliFitter1D cimport *
from String cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeFitter1D.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeFitter1D "OpenMS::IsotopeFitter1D":
        IsotopeFitter1D() nogil except +
        IsotopeFitter1D(IsotopeFitter1D) nogil except +
        # QualityType fit1d(RawDataArrayType &range, InterpolationModel *&model) nogil except +
        # Fitter1D * create() nogil except +
        String getProductName() nogil except +

