// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_MRMIONSERIES_H
#define OPENMS_ANALYSIS_OPENSWATH_MRMIONSERIES_H

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <boost/unordered_map.hpp>
#include <boost/assign.hpp>

// #define DEBUG_MRMIONSERIES

namespace OpenMS
{
  /**
    @brief Generate theoretical fragment ion series for use in MRMAssay and MRMDecoy

    Will generate theoretical fragment ionseries based on AASequence and parameters.
    Neutral losses are supported according to a model similar to the one in SpectraST.
    ReactionMonitoringTransition objects can be annotated with the corresponding CV
    terms.

  */
  class OPENMS_DLLAPI MRMIonSeries
  {
private:
    CVTermList annotationToCVTermList_(String annotation);

    void annotationToCV_(ReactionMonitoringTransition& tr);

public:
    //@{
    /// Constructor
    MRMIonSeries();

    /// Destructor
    ~MRMIonSeries();
    //@}

    typedef boost::unordered_map<String, double> IonSeries;

    /**
      @brief Selects ion from IonSeries according to annotation string

      @param ionseries the IonSeries from which to choose
      @param ionid the annotation string of the query fragment ion
      @value std::pair<String, double> the annotation and product m/z of
      the queried fragment ion

    */
    std::pair<String, double> getIon(IonSeries ionseries, String ionid);

    /**
      @brief Selects ion from IonSeries according to product m/z

      @param ionseries the IonSeries from which to choose
      @param ProductMz the product m/z of the queried fragment ion
      @param mz_threshold the m/z threshold for annotation of the fragment ion
      @value std::pair<String, double> the annotation and product m/z of
      the queried fragment ion

    */
    std::pair<String, double> annotateIon(IonSeries ionseries, double ProductMZ, double mz_threshold);

    /**
      @brief Annotates transition with CV terms

      @param tr the transition to annotate
      @param annotation the fragment ion annotation.

    */
    void annotateTransitionCV(ReactionMonitoringTransition& tr, String annotation);

    /**
      @brief Annotates transition

      @param tr the transition to annotate
      @param peptide the corresponding peptide
      @param mz_threshold the m/z threshold for annotation of the fragment ion
      @param enable_reannotation whether the original (e.g. SpectraST)
      annotation should be used or reannotation should be conducted
      @param fragment_types the fragment ion types for reannotation
      @param fragment_charges the fragment ion charges for reannotation
      @param enable_losses whether neutral losses should be considered
      @param round_decPow round product m/z values to decimal power (default: -4)

    */
    void annotateTransition(ReactionMonitoringTransition& tr, const TargetedExperiment::Peptide peptide, const double mz_threshold, bool enable_reannotation, std::vector<String> fragment_types, std::vector<size_t> fragment_charges, bool enable_losses, int round_decPow = -4);

    /**
      @brief Computed theoretical fragment ion series

      @param sequence the peptide amino acid sequence
      @param precursor_charge the charge of the peptide precursor
      @param fragment_types the fragment ion types for reannotation
      @param fragment_charges the fragment ion charges for reannotation
      @param enable_losses whether neutral losses should be considered
      @param round_decPow round product m/z values to decimal power (default: -4)
      @value IonSeries the theoretical fragment ion series
    */
    IonSeries getIonSeries(AASequence sequence, size_t precursor_charge, std::vector<String> fragment_types, std::vector<size_t> fragment_charges, bool enable_losses, int round_decPow = -4);
  };
}

#endif
