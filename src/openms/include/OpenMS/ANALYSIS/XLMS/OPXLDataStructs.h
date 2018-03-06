// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_XLMS_OPXLDATASTRUCTS
#define OPENMS_ANALYSIS_XLMS_OPXLDATASTRUCTS

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
//#include <numeric>

namespace OpenMS
{
  class OPENMS_DLLAPI OPXLDataStructs
  {

    public:

      /**
       * @brief The ProteinProteinCrossLinkType enum enumerates possible types of Protein-Protein cross-linking reaction results. Cross-link, Mono-link or Loop-link.
       */
      enum ProteinProteinCrossLinkType
      {
        CROSS = 0,
        MONO = 1,
        LOOP = 2,
        NUMBER_OF_CROSS_LINK_TYPES
      };

      /**
       * @brief The ProteinProteinCrossLink struct represents a cross-link between two peptides in OpenPepXL.

          This struct completely defines a cross-link linking two peptides. It contains the two peptides alpha and beta as AASequences,
          the positions that are linked in each peptide, the mass and name of the linker and terminal specificities to distinguish linkers attached to side chains
          of the terminal residues from linkers attached to the termini themselves.
          The beta peptide can be an empty AASequence object in case of a mono- or loop-link.
          Used to represent a theoretical candidate for generation of theoretical spectra and matching to experimental spectra.
       */
      struct ProteinProteinCrossLink
      {
        AASequence alpha; // longer peptide
        AASequence beta; // shorter peptide (empty for mono-link), tie bracker: mass then lexicographical
        std::pair<SignedSize, SignedSize> cross_link_position; // index in alpha, beta or between alpha, alpha in loop-links
        double cross_linker_mass;
        String cross_linker_name;
        ResidueModification::TermSpecificity term_spec_alpha;
        ResidueModification::TermSpecificity term_spec_beta;

        ProteinProteinCrossLinkType getType() const
        {
          if (!beta.empty()) return CROSS;

          if (cross_link_position.second == -1) return MONO;

          return LOOP;
        }

        bool operator==(const ProteinProteinCrossLink & other) const
        {
          return alpha == other.alpha &&
                      beta == other.beta &&
                      cross_link_position == other.cross_link_position &&
                      cross_linker_mass == other.cross_linker_mass &&
                      cross_linker_name == other.cross_linker_name &&
                      term_spec_alpha == other.term_spec_alpha &&
                      term_spec_beta == other.term_spec_beta;
        }
      };

      /**
       * @brief The CrossLinkSpectrumMatch struct represents a PSM between a ProteinProteinCrossLink and a spectrum in OpenPepXL.

          This struct contains a ProteinProteinCrossLink and indices to one or two spectra.
          It also contains the results of a match between the ProteinProteinCrossLink and these spectra as scores and peak annotations.
          Used as a temporary container to collect results efficiently, since only a few top matches will be kept for each experimental spectrum for output.
       */
      struct CrossLinkSpectrumMatch
      {
        /// structure of the cross-link
        ProteinProteinCrossLink cross_link;

        /// reference to pair of spectra
        Size scan_index_light;
        Size scan_index_heavy;

        /// final score
        double score;

        /// rank among the matches to the same spectrum
        Size rank;

        /// counts, scores and other data for xQuest-like output
        double pre_score;
        double percTIC;
        double wTIC;
        double int_sum;
        double match_odds;
        double log_occupancy;
        double log_occupancy_alpha;
        double log_occupancy_beta;
        double log_occupancy_full_spec;
        std::vector< double > xcorrx;
        double xcorrx_max;
        std::vector< double > xcorrc;
        double xcorrc_max;
        Size matched_common_alpha;
        Size matched_common_beta;
        Size matched_xlink_alpha;
        Size matched_xlink_beta;
        double HyperCommon;
        double HyperXlink;
        double HyperAlpha;
        double HyperBeta;
        double HyperBoth;
        double PScoreCommon;
        double PScoreXlink;
        double PScoreAlpha;
        double PScoreBeta;
        double PScoreBoth;

        std::vector<PeptideHit::PeakAnnotation> frag_annotations;

        Size peptide_id_index;

        bool operator<(const CrossLinkSpectrumMatch& other) const
        {
          return score < other.score;
        }

        bool operator==(const CrossLinkSpectrumMatch& other) const
        {
          return cross_link == other.cross_link &&
                     scan_index_light == other.scan_index_light &&
                     scan_index_heavy == other.scan_index_heavy &&
                     score == other.score &&
                     rank == other.rank &&
                     pre_score == other.pre_score &&
                     percTIC == other.percTIC &&
                     wTIC == other.wTIC &&
                     int_sum == other.int_sum &&
                     match_odds == other.match_odds &&
                     xcorrx == other.xcorrx &&
                     xcorrx_max == other.xcorrx_max &&
                     xcorrc == other.xcorrc &&
                     xcorrc_max == other.xcorrc_max &&
                     matched_common_alpha == other.matched_common_alpha &&
                     matched_common_beta == other.matched_common_beta &&
                     matched_xlink_alpha == other.matched_xlink_alpha &&
                     matched_xlink_beta == other.matched_xlink_beta &&
                     HyperCommon == other.HyperCommon &&
                     HyperXlink == other.HyperXlink &&
                     HyperAlpha == other.HyperAlpha &&
                     HyperBeta == other.HyperBeta &&
                     HyperBoth == other.HyperBoth &&
                     PScoreCommon == other.PScoreCommon &&
                     PScoreXlink == other.PScoreXlink &&
                     PScoreAlpha == other.PScoreAlpha &&
                     PScoreBeta == other.PScoreBeta &&
                     PScoreBoth == other.PScoreBoth;
        }
      };

      /**
       * @brief The XLPrecursor struct represents a cross-link candidate in the process of filtering candidates by precursor masses in OpenPepXL.

          Since the precursor mass does not change depending on the exact linked residues, one XLPrecursor can represent several ProteinProteinCrossLinks
          in the process of filtering by precursor mass. The precursor mass is the sum of the masses of the two peptides and the cross-linker.
          This struct also contains the indices of the two peptides in a vector, so that the two peptides can be identified again.
          This precursor mass is enumerated for all peptide pairs in the protein database given as input to OpenPepXL
          and is one of the major contributors to the memory usage of the tool because of the squared complexity of this enumeration.
          Therefore this should be kept as compact as possible.
       */
      struct XLPrecursor
      {
        float precursor_mass;
        unsigned int alpha_index;
        unsigned int beta_index;
      };

      // comparator for sorting XLPrecursor vectors and using upper_bound and lower_bound using only a precursor mass
      /**
       * @brief The XLPrecursorComparator is a comparator for XLPrecursors, that allows direct comparison of the XLPrecursor precursor mass with double numbers.

          This comparator can be used to sort a vector of XLPrecursors by precursor mass and to search for XLPrecursors with precursor masses within a double type mass range.
       */
      struct XLPrecursorComparator
      {
        bool operator() (const XLPrecursor& a, const XLPrecursor& b) const
        {
          return a.precursor_mass < b.precursor_mass;
        }
        bool operator() (const XLPrecursor& a, const double& b) const
        {
          return a.precursor_mass < b;
        }
        bool operator() (const double& a, const XLPrecursor& b) const
        {
          return a < b.precursor_mass;
        }
      };

      /**
       * @brief The PeptidePosition enum

          Used to record the positions of peptides in proteins after in silico digestion determine whether protein terminal modifications are possible on a peptide.
       */
      enum PeptidePosition
      {
        INTERNAL = 0,
        C_TERM = 1,
        N_TERM = 2
      };

      /**
       * @brief The AASeqWithMass struct represents a normal peptide with its precomputed mass.

          This struct stores information about a peptide as an AASequence and a PeptidePosition.
          It is used to enumerate pairs of peptides in OpenPepXL.
          Since the mass of every peptide is used many times, it is precomputed once and also stored in this struct.
          A vector of these structs is used to represent the digested protein database in OpenPepXL.
          An instance of this struct is created only once for each peptide in the digested database, so it does not contribute to memory usage
          as much as XLPrecursor does.
       */
      struct AASeqWithMass
      {
        double peptide_mass;
        AASequence peptide_seq;
        PeptidePosition position;
      };

      /**
       * @brief The AASeqWithMassComparator is a comparator for AASeqWithMass objects.

          This comparator allows to sort AASeqWithMass objects by the precomputed peptide mass and search for AASeqWithMass objects within a double type mass range.
       */
      struct AASeqWithMassComparator
      {
        bool operator() (const AASeqWithMass a, const AASeqWithMass b) const
        {
          return a.peptide_mass < b.peptide_mass;
        }
        bool operator() (const AASeqWithMass a, const double b) const
        {
          return a.peptide_mass < b;
        }
        bool operator() (const double a, const AASeqWithMass b) const
        {
          return a < b.peptide_mass;
        }
      };

      /**
       * @brief The PreprocessedPairSpectra struct represents the result of comparing a light and a heavy labeled spectra to each other.

          OpenPepXL can use labeled cross-linkers to denoise MS2 spectra. The PeakMaps contained in this struct represent the result of this
          denoising process for a whole mzML input file.
       */
      struct PreprocessedPairSpectra
      {

        MSExperiment spectra_common_peaks; // merge spectrum of common peaks (present in both spectra)
        MSExperiment spectra_xlink_peaks; // Xlink peaks in the light spectrum (common peaks between spectra_light_different and spectra heavy_to_light)
        MSExperiment spectra_all_peaks;

        // pre-initialize so we can simply std::swap the spectra (no synchronization in multi-threading context needed as we get no reallocation of the PeakMaps).
        PreprocessedPairSpectra(Size size)
        {
          for (Size i = 0; i != size; ++i)
          {
            spectra_common_peaks.addSpectrum(PeakSpectrum());
            spectra_xlink_peaks.addSpectrum(PeakSpectrum());
            spectra_all_peaks.addSpectrum(PeakSpectrum());
          }
        }
      };

  }; // class
} // namespace OpenMS


#endif
