from Types cimport *
from Matrix import *
from Feature cimport *
from FeatureMap cimport *
from String cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsotopeLabelingMDVs.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeLabelingMDVs(DefaultParamHandler) :
        
        IsotopeLabelingMDVs() nogil except +

        void isotopicCorrection(const Feature & normalized_feature, Feature & corrected_feature, 
          const Matrix[double] & correction_matrix, const DerivatizationAgent & correction_matrix_agent) nogil except +

        void isotopicCorrections(
          const FeatureMap & normalized_featureMap, FeatureMap & corrected_featureMap,
          const Matrix[double] & correction_matrix, const DerivatizationAgent & correction_matrix_agent) nogil except +

        void calculateIsotopicPurity(
          const Feature & normalized_feature, Feature & feature_with_isotopic_purity,
          const libcpp_vector[double] & experiment_data, const String & isotopic_purity_name) nogil except +

        void calculateIsotopicPurities(
          const FeatureMap & normalized_feature, FeatureMap & feature_with_isotopic_purity,
          const libcpp_vector[double] & experiment_data, const String & isotopic_purity_name) nogil except +

        void calculateMDVAccuracy(
          const Feature & normalized_feature, Feature & feature_with_accuracy_info,
          const libcpp_vector[double] & fragment_isotopomer_measured, const String & fragment_isotopomer_theoretical_formula) nogil except +

        void calculateMDVAccuracies(
          const FeatureMap & normalized_featureMap, FeatureMap & featureMap_with_accuracy_info,
          const libcpp_vector[double] & fragment_isotopomer_measured, const String & fragment_isotopomer_theoretical_formula) nogil except +

        void calculateMDV(
          const Feature & measured_feature, Feature & normalized_feature,
          const IsotopeLabelingMDVs_MassIntensityType & mass_intensity_type, const IsotopeLabelingMDVs_FeatureName & feature_name) nogil except +

        void calculateMDVs(
          const FeatureMap & measured_featureMap, FeatureMap & normalized_featureMap,
          const IsotopeLabelingMDVs_MassIntensityType & mass_intensity_type, const IsotopeLabelingMDVs_FeatureName & feature_name) nogil except +


cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsotopeLabelingMDVs.h>" namespace "OpenMS::IsotopeLabelingMDVs":

    cdef enum DerivatizationAgent "OpenMS::IsotopeLabelingMDVs::DerivatizationAgent":
        #wrap-attach:
        #    DerivatizationAgent
        NOT_SELECTED
        TBDMS
        SIZE_OF_DERIVATIZATIONAGENT

    cdef enum IsotopeLabelingMDVs_FeatureName "OpenMS::IsotopeLabelingMDVs::FeatureName":
        #wrap-attach:
        #    FeatureName
        INTENSITY
        PEAK_APEX_INT
        SIZE_OF_FEATURENAME

    cdef enum IsotopeLabelingMDVs_MassIntensityType "OpenMS::IsotopeLabelingMDVs::MassIntensityType":
        #wrap-attach:
        #    MassIntensityType
        NORM_MAX
        NORM_SUM
        SIZE_OF_MASSINTENSITYTYPE
