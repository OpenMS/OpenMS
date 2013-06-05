from InterpolationModel cimport *
from String cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>" namespace "OpenMS":
    
    cdef cppclass EmgModel(InterpolationModel):
        # wrap-inherits:
        #  InterpolationModel
        EmgModel() nogil except +
        EmgModel(EmgModel) nogil except +
        # BaseModel[ 1 ] * create() nogil except +
        String getProductName() nogil except +

        # inherited from parent class - no second definition necessary!
        # void setOffset(CoordinateType offset) # wrap-ignore
        # void setSamples() # wrap-ignore
        # CoordinateType getCenter() # wrap-ignore

