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

