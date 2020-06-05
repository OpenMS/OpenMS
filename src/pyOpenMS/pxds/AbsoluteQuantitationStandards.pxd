from Feature cimport *
from FeatureMap cimport *
from String cimport *
from Types cimport *

cdef extern from "<OpenMS/METADATA/AbsoluteQuantitationStandards.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitationStandards:

        AbsoluteQuantitationStandards() nogil except +
        AbsoluteQuantitationStandards(AbsoluteQuantitationStandards) nogil except +

        void getComponentFeatureConcentrations(
            libcpp_vector[AQS_runConcentration]& run_concentrations,
            libcpp_vector[FeatureMap]& feature_maps,
            const String& component_name,
            libcpp_vector[AQS_featureConcentration]& feature_concentrations
        ) nogil except +

        # void mapComponentsToConcentrations(libcpp_vector[ AQS_runConcentration ] run_concentrations,
        #                                    libcpp_vector[ FeatureMap ] feature_maps,
        #                                    libcpp_map[ String, libcpp_vector[ AQS_featureConcentration ]] & components_to_concentrations) nogil except +

cdef extern from "<OpenMS/METADATA/AbsoluteQuantitationStandards.h>" namespace "OpenMS::AbsoluteQuantitationStandards":

    cdef cppclass AQS_runConcentration "OpenMS::AbsoluteQuantitationStandards::runConcentration":

        AQS_runConcentration() nogil except +
        AQS_runConcentration(AQS_runConcentration) nogil except +

        String sample_name
        String component_name
        String IS_component_name
        double actual_concentration
        double IS_actual_concentration
        String concentration_units
        double dilution_factor

    cdef cppclass AQS_featureConcentration "OpenMS::AbsoluteQuantitationStandards::featureConcentration":

        AQS_featureConcentration() nogil except +
        AQS_featureConcentration(AQS_featureConcentration) nogil except +

        Feature feature
        Feature IS_feature
        double actual_concentration
        double IS_actual_concentration
        String concentration_units
        double dilution_factor

