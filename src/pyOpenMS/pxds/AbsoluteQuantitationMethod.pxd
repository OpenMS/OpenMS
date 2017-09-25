from Types cimport *
from Param cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitationMethod:

        AbsoluteQuantitationMethod()  nogil except +
        AbsoluteQuantitationMethod(AbsoluteQuantitationMethod)  nogil except + #wrap-ignore

        void setLOD(double llod, double ulod) nogil except +
        # void getLOD(double llod, double ulod) nogil except +
        void setLOQ(double lloq, double uloq) nogil except +
        # void getLOQ(double lloq, double uloq) nogil except +
        void checkLOD(double value) nogil except +
        void checkLOQ(double value) nogil except +
        void setComponentISFeatureNames(String component_name, String IS_name, String feature_name) nogil except +
        void getComponentISFeatureNames(String component_name, String IS_name, String feature_name) nogil except +
        void setConcentrationUnits(String concentration_units) nogil except +
        void getConcentrationUnits(String concentration_units) nogil except +
        void setTransformationModel(String transformation_model, Param & transformation_model_param) nogil except +
        void getTransformationModel(String transformation_model, Param & transformation_model_param) nogil except +
        void setActualConcentration(double actual_concentration) nogil except +
        # void getActualConcentration(double actual_concentration) nogil except +
        void setStatistics(int n_points, double correlation_coefficient) nogil except +
        # void getStatistics(int n_points, double correlation_coefficient) nogil except +
        Param fitTransformationModel(String transformation_model, libcpp_vector[ libcpp_pair[double, double] ]& data, Param & transformation_model_params) nogil except +
        double evaluateTransformationModel(String transformation_model, double datum, Param & transformation_model_params) nogil except +
