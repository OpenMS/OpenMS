from InterpolationModel cimport *
from String cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>" namespace "OpenMS":
    
    cdef cppclass EmgModel(InterpolationModel):
        # wrap-inherits:
        #  InterpolationModel
        EmgModel() except + nogil  # wrap-doc:Exponentially modified gaussian distribution model for elution profiles
        EmgModel(EmgModel &) except + nogil 
        # BaseModel[ 1 ] * create() except + nogil 
        String getProductName() except + nogil  # wrap-doc:Name of the model 

        # inherited from parent class - no second definition necessary!
        # void setOffset(CoordinateType offset) # wrap-ignore
        # void setSamples() # wrap-ignore
        # CoordinateType getCenter() # wrap-ignore
