from Types cimport *
from Param cimport *
from Feature import *

cdef extern from "<OpenMS/METADATA/AbsoluteQuantitationStandards.h>" namespace "OpenMS":    

    cdef cppclass AQS_runConcentration "OpenMS::AbsoluteQuantitationStandards::runConcentration":
        AQS_runConcentration()
        AQS_runConcentration(AQS_runConcentration &) # no-wrap

        String run_id
        String component_id
        String IS_component_id
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
