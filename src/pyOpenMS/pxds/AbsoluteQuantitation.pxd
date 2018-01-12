from Types cimport *
from Param cimport *
from AbsoluteQuantitationMethod cimport *
from AbsoluteQuantitationStandards cimport *
from Feature cimport *
from FeatureMap cimport *
from String cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitation(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler

        AbsoluteQuantitation() nogil except +
        AbsoluteQuantitation(AbsoluteQuantitation)  nogil except + #wrap-ignore

        void setQuantMethods(libcpp_vector[ AbsoluteQuantitationMethod ]& quant_methods) nogil except +
        libcpp_vector[ AbsoluteQuantitationMethod ] getQuantMethods() nogil except +
        double calculateRatio(Feature & component_1, Feature & component_2, String & feature_name) nogil except +
        # double calculateBias(double actual_concentration, double calculated_concentration) nogil except +
        double applyCalibration(Feature & component, Feature & IS_component, String & feature_name, 
                                String & transformation_model, Param & transformation_model_params) nogil except +
        void quantifyComponents(FeatureMap& unknowns) nogil except +

        void optimizeCalibrationCurveIterative(
            libcpp_vector[ AQS_featureConcentration ] & component_concentrations,
            String & feature_name,
            String & transformation_model,
            Param & transformation_model_params,
            Param & optimized_params) nogil except +

        void optimizeCalibrationCurves(
           libcpp_map[ String, libcpp_vector[ AQS_featureConcentration ]] & components_concentrations) nogil except + # wrap-ignore


