# from MaxLikeliFitter1D cimport *
from String cimport *

cdef extern from "<OpenMS/FEATUREFINDER/IsotopeFitter1D.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeFitter1D "OpenMS::IsotopeFitter1D":

        IsotopeFitter1D() except + nogil  # wrap-doc:Isotope distribution fitter (1-dim.) approximated using linear interpolation
        IsotopeFitter1D(IsotopeFitter1D &) except + nogil 

        # QualityType fit1d(RawDataArrayType &range, InterpolationModel *&model) except + nogil 
        # Fitter1D * create() except + nogil 

