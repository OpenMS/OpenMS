// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <unordered_map>

// #define DEBUG_MRMIONSERIES

namespace OpenMS
{

  /**
    @brief Generate theoretical fragment ion series for use in MRMAssay and MRMDecoy

    Will generate theoretical fragment ionseries based on AASequence and parameters.
    Neutral losses are supported according to a model similar to the one in SpectraST.
    ReactionMonitoringTransition objects can be annotated with the corresponding CV
    terms.

    MRMIonSeries uses internally an annotation format that is compatible with SpectraST,
    which is derived from Roepstoff and Fohlman (1984, PMID: 6525415). An annotation tag
    starts with ion type (a, b, c, x, y, z) followed by the ordinal (number of amino acids).
    The caret symbol follows "^" with a positive integer indicating the fragment ion charge.
    If no caret symbol is present, a charge of 1 is assumed. In case of neutral loss, a
    negative symbol "-" followed by the integer mass (e.g. 17 for ammonia) OR the molecular
    composition, compatible with EmpiricalFormula (e.g. N1H3 for ammonia) is allowed.

    Valid examples: y3, y3^1, y3^1-18, y3^1-H2O, y3-H2O

    Limitations: Special SpectraST multi-assignments, immonium, precursors are not supported.

  */
  class OPENMS_DLLAPI MRMIonSeries
  {
private:
    TargetedExperiment::Interpretation annotationToCVTermList_(const String& annotation);

    void annotationToCV_(ReactionMonitoringTransition& tr);

public:
    //@{
    /// Constructor
    MRMIonSeries();

    /// Destructor
    ~MRMIonSeries();
    //@}

    typedef std::unordered_map<String, double> IonSeries; ///< An MRM ion series which maps: "ion_type" -> "fragment m/z"

    /**
      @brief Selects ion from IonSeries according to annotation string

      @param[in,out] ionseries the IonSeries from which to choose
      @param[in] ionid the annotation string of the query fragment ion
      @return std::pair<String, double> the annotation and product m/z of
      the queried fragment ion

    */
    std::pair<String, double> getIon(IonSeries& ionseries, const String& ionid);

    /**
      @brief Selects ion from IonSeries according to product m/z

      @param[in] ionseries the IonSeries from which to choose
      @param[in] product_mz the product m/z of the queried fragment ion
      @param[in] mz_threshold the m/z threshold for annotation of the fragment ion
      @return std::pair<String, double> the annotation and product m/z of
      the queried fragment ion

    */
    std::pair<String, double> annotateIon(const IonSeries& ionseries, const double product_mz, const double mz_threshold);

    /**
      @brief Annotates transition with CV terms

      @param[in,out] tr the transition to annotate
      @param[in] annotation the fragment ion annotation.

    */
    void annotateTransitionCV(ReactionMonitoringTransition& tr, const String& annotation);

    /**
      @brief Annotates transition

      @param[in,out] tr the transition to annotate
      @param[in] peptide the corresponding peptide
      @param[in] precursor_mz_threshold the m/z threshold for annotation of the precursor ion
      @param[in] product_mz_threshold the m/z threshold for annotation of the fragment ion
      @param[in] enable_reannotation whether the original (e.g. SpectraST)
      annotation should be used or reannotation should be conducted
      @param[in] fragment_types the fragment ion types for reannotation
      @param[in] fragment_charges the fragment ion charges for reannotation
      @param[in] enable_specific_losses whether specific neutral losses should be considered
      @param[in] enable_unspecific_losses whether unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) should be considered
      @param[in] round_decPow round precursor and product m/z values to decimal power (default: -4)

    */
    void annotateTransition(ReactionMonitoringTransition& tr,
                            const TargetedExperiment::Peptide& peptide,
                            const double precursor_mz_threshold,
                            const double product_mz_threshold,
                            const bool enable_reannotation,
                            const std::vector<String>& fragment_types,
                            const std::vector<size_t>& fragment_charges,
                            const bool enable_specific_losses,
                            const bool enable_unspecific_losses,
                            const int round_decPow = -4);

    /**
      @brief Computed theoretical fragment ion series

      @param[in] sequence the peptide amino acid sequence
      @param[in] precursor_charge the charge of the peptide precursor
      @param[in] fragment_types the fragment ion types for reannotation
      @param[in] fragment_charges the fragment ion charges for reannotation
      @param[in] enable_specific_losses whether specific neutral losses should be considered
      @param[in] enable_unspecific_losses whether unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) should be considered
      @param[in] round_decPow round product m/z values to decimal power (default: -4)
      @return IonSeries the theoretical fragment ion series
    */
    IonSeries getIonSeries(const AASequence& sequence,
                           size_t precursor_charge,
                           const std::vector<String>& fragment_types,
                           const std::vector<size_t>& fragment_charges,
                           const bool enable_specific_losses,
                           const bool enable_unspecific_losses,
                           const int round_decPow = -4);
  };
}

