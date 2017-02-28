from Types cimport *
from MZTrafoModel cimport *
from RANSAC cimport *
from String cimport *
from CalibrationData cimport *

cdef extern from "<OpenMS/FILTERING/CALIBRATION/MZTrafoModel.h>" namespace "OpenMS":

    cdef cppclass MZTrafoModel:
        MZTrafoModel()  nogil except +
        MZTrafoModel(MZTrafoModel &) nogil except +
        MZTrafoModel(bool) nogil except +
        bool isTrained() nogil except +
        double getRT() nogil except +
        double predict(double) nogil except +
        bool train(CalibrationData, MZTrafoModel_MODELTYPE, bool, double, double) nogil except +
        bool train(libcpp_vector[double], libcpp_vector[double], libcpp_vector[double], MZTrafoModel_MODELTYPE, bool) nogil except +
        void getCoefficients(double& intercept, double& slope, double& power) nogil except +
        void setCoefficients(MZTrafoModel) nogil except +
        void setCoefficients(double, double, double) nogil except +
        String toString() nogil except +
        
cdef extern from "<OpenMS/FILTERING/CALIBRATION/MZTrafoModel.h>" namespace "OpenMS::MZTrafoModel":

    cdef enum MZTrafoModel_MODELTYPE "OpenMS::MZTrafoModel::MODELTYPE":
        LINEAR
        LINEAR_WEIGHTED
        QUADRATIC
        QUADRATIC_WEIGHTED
        SIZE_OF_MODELTYPE
    
    # static members
    # libcpp_string names_of_modeltype[] nogil except +
    MZTrafoModel_MODELTYPE nameToEnum(libcpp_string name) nogil except + # wrap-attach:MZTrafoModel
    libcpp_string enumToName(MZTrafoModel_MODELTYPE mt) nogil except + # wrap-attach:MZTrafoModel
    void setRANSACParams(RANSACParam p) nogil except + # wrap-attach:MZTrafoModel
    void setCoefficientLimits(double offset, double scale, double power) nogil except + # wrap-attach:MZTrafoModel
    bool isValidModel(MZTrafoModel& trafo) nogil except + # wrap-attach:MZTrafoModel
    Size findNearest(libcpp_vector[MZTrafoModel]& tms, double rt) nogil except + # wrap-attach:MZTrafoModel
    
