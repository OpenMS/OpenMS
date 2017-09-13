from Types cimport *
from Param cimport *
from AbsoluteQuantitationMethod cimport *
from AbsoluteQuantitationStandards cimport *
from Feature cimport *
from FeatureMap cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitation:
        # wrap-ignore
        # no-pxd-import

        AbsoluteQuantitation() nogil except +

        void setQuantMethods(libcpp_vector[ AbsoluteQuantitationMethod ]& quant_methods) nogil except +
        double calculateRatio(Feature & component_1,Feature & component_2, string feature_name) nogil except +
        double calculateBias(double actual_concentration, double calculated_concentration) nogil except +
        double applyCalibration(Feature & component,Feature & IS_component, string feature_name, string transformation_model, Param & transformation_model_params) nogil except +
        void quantifyComponents(libcpp_vector[ FeatureMap ]& unknowns) nogil except +

