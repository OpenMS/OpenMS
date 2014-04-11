from InterpolationModel cimport *
from IsotopeDistribution cimport *
from EmpiricalFormula cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeModel "OpenMS::IsotopeModel":
        IsotopeModel() nogil except +
        IsotopeModel(IsotopeModel) nogil except +
        UInt getCharge() nogil except +
        void setOffset(double offset) nogil except +
        double getOffset() nogil except +
        EmpiricalFormula getFormula() nogil except +
        void setSamples(EmpiricalFormula &formula) nogil except +
        double getCenter() nogil except +
        IsotopeDistribution  getIsotopeDistribution() nogil except +
        # BaseModel[ 1 ] * create() nogil except +
        String getProductName() nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>" namespace "OpenMS::IsotopeModel":
    cdef enum Averagines "OpenMS::IsotopeModel::Averagines":
        #wrap-attach:
        #    IsotopeModel
        C
        H
        N
        O
        S
        AVERAGINE_NUM

