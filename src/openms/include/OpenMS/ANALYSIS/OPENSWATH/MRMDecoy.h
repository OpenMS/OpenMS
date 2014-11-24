// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_MRMDECOY_H
#define OPENMS_ANALYSIS_OPENSWATH_MRMDECOY_H

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/assign.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/unordered_map.hpp>

#include <string>
#include <vector>
#include <map>
#include <utility> // for pair

// #define DEBUG_MRMDECOY

namespace OpenMS
{
/**
  @brief This class generates a TargetedExperiment object with decoys based on a TargetedExperiment object
 
  There are multiple methods to create the decoy transitions, the simplest ones
  are reverse and pseudo-reverse which reverse the sequence either completely or
  leaving the last (tryptic) AA untouched respectively.
 
  Another decoy generation method is "shuffle" which uses an algorithm similar
  to the one described in Lam, Henry, et al. (2010). "Artificial decoy spectral
  libraries for false discovery rate estimation in spectral library searching in
  proteomics".  Journal of Proteome Research 9, 605-610. It shuffles the amino
  acid sequence and shuffles the fragment ion intensities accordingly, however
  for this to work the fragment ions need to be matched to annotated before.


  First, the algorithm goes through all peptides and applies the decoy method to
  the target peptide sequence (pseudo-reverse, reverse or shuffle) in order to
  produce the decoy sequence. Then, for each peptide, the fragment ions in the
  target library are matched to their most likely origin (e.g. the ions are
  annotated with their ion series (a,b,y) and the fragment number and optionally
  a neutral loss (10 different neutral losses are currently implemented)). For
  each fragment ion from the target peptide, an equivalent ion is created for the
  decoy peptide with the same intensity (e.g. if the target peptide sequence has
  a b5 ion with a normalized intensity of 200, an equivalent b5 ion for the
  decoy sequence is created and assigned the intensity 200).
  Optionally, the m/z values are corrected to reflect the theoretical value rather
  than the experimental value in the library.

 */

  class OPENMS_DLLAPI MRMDecoy :
    public ProgressLogger
  {
public:
    MRMDecoy() {} // empty, no members

    /**
      @brief Generate decoys from a TargetedExperiment

      Will generate decoy peptides for each target peptide provided in exp and
      write them into the decoy experiment.

      Valid methods: shuffle, reverse, pseudo-reverse

      If theoretical is true, the target transitions will be returned but their
      masses will be adjusted to match the theoretical value of the fragment ion
      that is the most likely explanation for the product.

      mz_threshold is used for the matching of theoretical ion series to the observed one

    */
    void generateDecoys(OpenMS::TargetedExperiment& exp,
                        OpenMS::TargetedExperiment& dec, String method, String decoy_tag,
                        double identity_threshold, int max_attempts, double mz_threshold, 
                        bool theoretical, double mz_shift, bool exclude_similar, 
                        double similarity_threshold, bool remove_CNterm_mods, 
                        double precursor_mass_shift, bool enable_losses, bool remove_unannotated);
    /**
      @brief Remove transitions s.t. all peptides have a defined set of transitions.

      All transitions of a peptide above max_transitions get deleted, all
      peptides with less than min_transitions also get deleted.

    */
    void restrictTransitions(OpenMS::TargetedExperiment& exp, int min_transitions,
                             int max_transitions);

    typedef std::vector<OpenMS::TargetedExperiment::Protein> ProteinVectorType;
    typedef std::vector<OpenMS::TargetedExperiment::Peptide> PeptideVectorType;
    typedef std::vector<OpenMS::ReactionMonitoringTransition> TransitionVectorType;

    typedef boost::unordered_map<String, boost::unordered_map<String, double> > IonSeries;
    typedef boost::unordered_map<String, IonSeries> IonSeriesMapType;

    typedef std::map<String, std::vector<const ReactionMonitoringTransition*> > PeptideTransitionMapType;

    /**
      @brief Selects a decoy ion from a set of ions.
    */
    std::pair<String, double> getDecoyIon(String ionid,
                                              boost::unordered_map<String, boost::unordered_map<String, double> >& decoy_ionseries);

    /**
      @brief Selects a target ion from a set of ions.
    */
    std::pair<String, double> getTargetIon(double ProductMZ, double mz_threshold,
                                           boost::unordered_map<String, boost::unordered_map<String, double> > target_ionseries,
                                           bool enable_losses);
    /**
      @brief Generate all ion series for an input AASequence

      Currently generated are:

      bionseries, bionseries_loss,
      yionseries, yionseries_loss,
      aionseries

      for each of these, the following neutral losses are calculated:
        -17, -18, -34, -35, -36, -44, -45, -46, -64, -98.

      FEATURE (george): a more generic mechanism to specify which series and losses should be
      generated. possible integration with TheoreticalSpectrumGenerator?
    */
    boost::unordered_map<String, boost::unordered_map<String, double> > getIonSeries(
      AASequence sequence, int precursor_charge);

    /**
      @brief Find all tryptic sites in a sequence
    */
    std::vector<std::pair<std::string::size_type, std::string> > find_all_tryptic(
      std::string sequence);

    /**
      @brief Compute relative identity (relative number of matches of amino acids at the same position) between two sequences
    */
    float AASequenceIdentity(const String& sequence, const String& decoy);

    /**
      @brief Check if a peptide has C or N terminal modifications 
    */
    bool has_CNterminal_mods(const OpenMS::TargetedExperiment::Peptide & peptide);

    /**
      @brief Correct the masses according to theoretically computed masses
    */
    void correctMasses(OpenMS::TargetedExperiment& exp, double mz_threshold, bool enable_losses);

    /**
      @brief Shuffle a peptide (with its modifications) sequence

      This function will shuffle the given peptide sequences and its
      modifications such that the resulting relative sequence identity is below
      identity_threshold.
    */
    OpenMS::TargetedExperiment::Peptide shufflePeptide(
      OpenMS::TargetedExperiment::Peptide peptide, double identity_threshold, int seed = -1,
      int max_attempts = 10);

    /**
      @brief Pseudo-reverse a peptide sequence (with its modifications)

      Pseudo reverses a peptide sequence, leaving the last AA constant
    */
    OpenMS::TargetedExperiment::Peptide pseudoreversePeptide(
      OpenMS::TargetedExperiment::Peptide peptide);

    /**
      @brief Reverse a peptide sequence (with its modifications)
    */
    OpenMS::TargetedExperiment::Peptide reversePeptide(
      OpenMS::TargetedExperiment::Peptide peptide);
  };
}

#endif
