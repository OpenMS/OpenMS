from Types cimport *
from Matrix import *
from Feature cimport *
from FeatureMap cimport *
from String cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsotopeLabelingMDVs.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeLabelingMDVs() :
        
        IsotopeLabelingMDVs() nogil except +

        void isotopicCorrection(const FeatureMap & normalized_featureMap, FeatureMap & corrected_featureMap, 
          const MatrixDouble & correction_matrix, const String correction_matrix_agent) nogil except +

        void isotopicCorrections(
          const FeatureMap & normalized_featureMap, FeatureMap & corrected_featureMap,
          const MatrixDouble & correction_matrix, const String correction_matrix_agent) nogil except +

        void calculateIsotopicPurity(
          const Feature & normalized_feature, Feature & feature_with_isotopic_purity,
          libcpp_vector[double] & experiment_data, String & isotopic_purity_name) nogil except +

        void calculateIsotopicPurities(
          const FeatureMap & normalized_featureMap, FeatureMap & featureMap_with_isotopic_purity,
          libcpp_vector[double] & experiment_data, String & isotopic_purity_name) nogil except +

        void calculateMDVAccuracy(
          const Feature & normalized_feature, Feature & feature_with_accuracy_info,
          const libcpp_vector[double] & fragment_isotopomer_measured, const libcpp_vector[double] & fragment_isotopomer_theoretical) nogil except +

        void calculateMDVAccuracies(
          const FeatureMap & normalized_featureMap, FeatureMap & featureMap_with_accuracy_info,
          const libcpp_vector[double] & fragment_isotopomer_measured, const libcpp_vector[double] & fragment_isotopomer_theoretical) nogil except +

        void calculateMDV(
          const Feature & measured_feature, Feature & normalized_feature,
          const String & mass_intensity_type, const String & feature_name) nogil except +

        void calculateMDVs(
          const FeatureMap & measured_featureMap, FeatureMap & normalized_featureMap,
          const String & mass_intensity_type, const String & feature_name) nogil except +
