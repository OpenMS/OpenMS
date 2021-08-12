# from MaxLikeliFitter1D cimport *
from String cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeFitter1D.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeFitter1D "OpenMS::IsotopeFitter1D":

        IsotopeFitter1D() nogil except + # wrap-doc:Isotope distribution fitter (1-dim.) approximated using linear interpolation
        IsotopeFitter1D(IsotopeFitter1D &) nogil except +

        # QualityType fit1d(RawDataArrayType &range, InterpolationModel *&model) nogil except +
        # Fitter1D * create() nogil except +
        String getProductName() nogil except + # wrap-doc:Name of the model (needed by Factory)

