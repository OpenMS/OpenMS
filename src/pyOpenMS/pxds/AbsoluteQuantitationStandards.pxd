from Types cimport *
from Param cimport *
from Feature cimport *

cdef extern from "<OpenMS/METADATA/AbsoluteQuantitationStandards.h>" namespace "OpenMS":    

    cdef cppclass AQS_runConcentration "OpenMS::AbsoluteQuantitationStandards::runConcentration":
        AQS_runConcentration()
        AQS_runConcentration(AQS_runConcentration &) # no-wrap

        String sample_name
        String component_name
        String IS_component_name
        double actual_concentration
        double IS_actual_concentration
        String concentration_units
        double dilution_factor

    cdef cppclass AQS_featureConcentration "OpenMS::AbsoluteQuantitationStandards::featureConcentration":
        AQS_featureConcentration()
        AQS_featureConcentration(AQS_featureConcentration &) # no-wrap

        Feature feature
        Feature IS_feature
        double actual_concentration
        double IS_actual_concentration
        String concentration_units
        double dilution_factor

    cdef cppclass AbsoluteQuantitationStandards:

        AbsoluteQuantitationStandards() nogil except +
        AbsoluteQuantitationStandards(AbsoluteQuantitationStandards)  nogil except + #wrap-ignore
