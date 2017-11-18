from Types cimport *
from Param cimport *
from AbsoluteQuantitationMethod cimport *
from AbsoluteQuantitationStandards cimport *
from Feature cimport *
from FeatureMap cimport *
from String cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitation:

        AbsoluteQuantitation() nogil except +
        AbsoluteQuantitation(AbsoluteQuantitation)  nogil except + #wrap-ignore

        void setQuantMethods(libcpp_vector[ AbsoluteQuantitationMethod ]& quant_methods) nogil except +
        double calculateRatio(Feature & component_1, Feature & component_2, String & feature_name) nogil except +
        # double calculateBias(double actual_concentration, double calculated_concentration) nogil except +
        double applyCalibration(Feature & component, Feature & IS_component, String & feature_name, 
                                String & transformation_model, Param & transformation_model_params) nogil except +
        void quantifyComponents(FeatureMap& unknowns) nogil except +

