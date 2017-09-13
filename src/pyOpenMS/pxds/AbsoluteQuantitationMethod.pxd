from Types cimport *
from Param cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitationMethod:
        # wrap-ignore
        # no-pxd-import

        AbsoluteQuantitationMethod()  nogil except +

        void setLOD(double llod, double ulod) nogil except +
        void getLOD(double llod, double ulod) nogil except +
        void setLOQ(double lloq, double uloq) nogil except +
        void getLOQ(double lloq, double uloq) nogil except +
        void checkLOD(double value) nogil except +
        void checkLOQ(double value) nogil except +
        void setComponentISFeatureNames(string component_name, string IS_name, string, feature_name) nogil except +
        void getComponentISFeatureNames(string component_name, string IS_name, string, feature_name) nogil except +
        void setConcentrationUnits(string concentration_units) nogil except +
        void getConcentrationUnits(string concentration_units) nogil except +
        void setTransformationModel(string transformation_model, Param & transformation_model_param) nogil except +
        void getTransformationModel(string transformation_model, Param & transformation_model_param) nogil except +
        void setActualConcentration(double actual_concentration) nogil except +
        void getActualConcentration(double actual_concentration) nogil except +
        void setStatistics(int n_points, double correlation_coefficient) nogil except +
        void getStatistics(int n_points, double correlation_coefficient) nogil except +
        Param fitTransformationModel(string transformation_model, libcpp_vector[ libcpp_pair[double, double] ]& data, Param & transformation_model_params) nogil except +
        double evaluateTransformationModel(string transformation_model, double datum, Param & transformation_model_params) nogil except +