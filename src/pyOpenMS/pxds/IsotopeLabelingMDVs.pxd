from Types cimport *
from Matrix cimport *
from libcpp.map cimport map as libcpp_map
from DoubleList cimport *
from String cimport *
from Feature cimport *
from FeatureMap cimport *
from String cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsotopeLabelingMDVs.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeLabelingMDVs(DefaultParamHandler) :
        
        IsotopeLabelingMDVs() nogil except +
        IsotopeLabelingMDVs(IsotopeLabelingMDVs &) nogil except + # compiler

        void isotopicCorrection(const Feature & normalized_feature, Feature & corrected_feature, 
          Matrix[double] & correction_matrix, const DerivatizationAgent & correction_matrix_agent) nogil except +
            # wrap-doc:
              #   This function performs an isotopic correction to account for unlabeled abundances coming from
              #   the derivatization agent (e.g., tBDMS) using correction matrix method and is calculated as follows:
              #   -----
              #   :param  normalized_feature: Feature with normalized values for each component and unlabeled chemical formula for each component group
              #   :param  correction_matrix: Square matrix holding correction factors derived either experimentally or theoretically which describe how spectral peaks of naturally abundant 13C contribute to spectral peaks that overlap (or convolve) the spectral peaks of the corrected MDV of the derivatization agent
              #   :param  correction_matrix_agent: Name of the derivatization agent, the internally stored correction matrix if the name of the agent is supplied, only "TBDMS" is supported for now
              #   :returns: corrected_feature: Feature with corrected values for each component

        void isotopicCorrections(
          const FeatureMap & normalized_featureMap, FeatureMap & corrected_featureMap,
          Matrix[double] & correction_matrix, const DerivatizationAgent & correction_matrix_agent) nogil except +
            # wrap-doc:
              #   This function performs an isotopic correction to account for unlabeled abundances coming from
              #   the derivatization agent (e.g., tBDMS) using correction matrix method and is calculated as follows:
              #   -----
              #   :param  normalized_featuremap: FeatureMap with normalized values for each component and unlabeled chemical formula for each component group
              #   :param  correction_matrix: Square matrix holding correction factors derived either experimentally or theoretically which describe how spectral peaks of naturally abundant 13C contribute to spectral peaks that overlap (or convolve) the spectral peaks of the corrected MDV of the derivatization agent
              #   :param  correction_matrix_agent: Name of the derivatization agent, the internally stored correction matrix if the name of the agent is supplied, only "TBDMS" is supported for now
              #   :returns corrected_featuremap: FeatureMap with corrected values for each component

        void calculateIsotopicPurity(
          const Feature & normalized_feature,
          const libcpp_vector[double] & experiment_data, const String & isotopic_purity_name) nogil except +
            # wrap-doc:
              #   This function calculates the isotopic purity of the MDV using the following formula:
              #   isotopic purity of tracer (atom % 13C) = n / [n + (M + n-1)/(M + n)],
              #   where n in M+n is represented as the index of the result
              #   The formula is extracted from "High-resolution 13C metabolic flux analysis",
              #   Long et al, doi:10.1038/s41596-019-0204-0
              #   -----
              #   :param normalized_feature: Feature with normalized values for each component and the number of heavy labeled e.g., carbons. Out is a Feature with the calculated isotopic purity for the component group
              #   :param experiment_data: Vector of experiment data in percent
              #   :param isotopic_purity_name: Name of the isotopic purity tracer to be saved as a meta value

        # void calculateIsotopicPurities(
        #  const FeatureMap & normalized_feature,
        #  const libcpp_vector[ DoubleList ] & experiment_data, const libcpp_vector[String] & isotopic_purity_name) nogil except +

        void calculateMDVAccuracy(
          const Feature & normalized_feature,
          const String & feature_name, const String & fragment_isotopomer_theoretical_formula) nogil except +
            # wrap-doc:
              #   This function calculates the accuracy of the MDV as compared to the theoretical MDV (only for 12C quality control experiments)
              #   using average deviation to the mean. The result is mapped to the meta value "average_accuracy" in the updated feature
              #   -----
              #   :param normalized_feature: Feature with normalized values for each component and the chemical formula of the component group. Out is a Feature with the component group accuracy and accuracy for the error for each component
              #   :param fragment_isotopomer_measured: Measured scan values
              #   :param fragment_isotopomer_theoretical_formula: Empirical formula from which the theoretical values will be generated

        void calculateMDVAccuracies(
          const FeatureMap & normalized_featureMap,
          const String & feature_name, const libcpp_map[ libcpp_string, libcpp_string ] & fragment_isotopomer_theoretical_formulas) nogil except +
            # wrap-doc:
              #   This function calculates the accuracy of the MDV as compared to the theoretical MDV (only for 12C quality control experiments)
              #   using average deviation to the mean
              #   -----
              #   param normalized_featuremap: FeatureMap with normalized values for each component and the chemical formula of the component group. Out is a FeatureMap with the component group accuracy and accuracy for the error for each component
              #   param fragment_isotopomer_measured: Measured scan values
              #   param fragment_isotopomer_theoretical_formula: A map of ProteinName/peptideRef to Empirical formula from which the theoretical values will be generated

        void calculateMDV(
          const Feature & measured_feature, Feature & normalized_feature,
          const MassIntensityType & mass_intensity_type, const String & feature_name) nogil except +

        void calculateMDVs(
          const FeatureMap & measured_featureMap, FeatureMap & normalized_featureMap,
          const MassIntensityType & mass_intensity_type, const String & feature_name) nogil except +


cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsotopeLabelingMDVs.h>" namespace "OpenMS::IsotopeLabelingMDVs":

    cdef enum DerivatizationAgent "OpenMS::IsotopeLabelingMDVs::DerivatizationAgent":
        #wrap-attach:
        #    DerivatizationAgent
        NOT_SELECTED,
        TBDMS,
        SIZE_OF_DERIVATIZATIONAGENT

    cdef enum MassIntensityType "OpenMS::IsotopeLabelingMDVs::MassIntensityType":
        #wrap-attach:
        #    MassIntensityType
        NORM_MAX,
        NORM_SUM,
        SIZE_OF_MASSINTENSITYTYPE
