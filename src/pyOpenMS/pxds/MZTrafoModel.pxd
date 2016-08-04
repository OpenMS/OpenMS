from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from MZTrafoModel cimport *

cdef extern from "<OpenMS/FILTERING/CALIBRATION/MZTrafoModel.h>" namespace "OpenMS":

    cdef cppclass MZTrafoModel:
        MZTrafoModel()  nogil except +
        MZTrafoModel(MZTrafoModel &) nogil except +
        MZTrafoModel(bool) nogil except +
        bool isTrained() nogil except +
        double getRT() nogil except +
        double predict(double) nogil except +
        bool train(CalibrationData, MODELTYPE, bool, double, double) nogil except +
        bool train(libcpp_vector[double], libcpp_vector[double], libcpp_vector[double], MODELTYPE, bool) nogil except +
        # return by ref:
        #void getCoefficients(double& intercept, double& slope, double& power) nogil except +
        void setCoefficients(MZTrafoModel) nogil except +
        void setCoefficients(double, double, double) nogil except +
        String toString() nogil except +
        
cdef extern from "<OpenMS/FILTERING/CALIBRATION/MZTrafoModel.h>" namespace "OpenMS::MZTrafoModel":

    cdef enum MODELTYPE "OpenMS::MZTrafoModel::MODELTYPE":
        LINEAR
        LINEAR_WEIGHTED
        QUADRATIC
        QUADRATIC_WEIGHTED
        SIZE_OF_MODELTYPE
    
    # static members
    const std::string names_of_modeltype[] nogil except +
    MODELTYPE nameToEnum(const std::string& name) nogil except +
    const std::string& enumToName(MODELTYPE mt) nogil except + 
    void setRANSACParams(const Math::RANSACParam& p) nogil except +
    void setCoefficientLimits(double offset, double scale, double power) nogil except +
    bool isValidModel(const MZTrafoModel& trafo) nogil except +
    Size findNearest(const std::vector<MZTrafoModel>& tms, double rt) nogil except +

    