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

#ifndef OPENMS_ANALYSIS_XLMS_OPENPROXLUTILS
#define OPENMS_ANALYSIS_XLMS_OPENPROXLUTILS

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <numeric>

namespace OpenMS
{

  struct ProteinProteinCrossLink
    {
      /** Enums
      */
      //@{
      /** @brief type of Protein-Protein cross-link
      */
      enum ProteinProteinCrossLinkType
      {
        CROSS = 0,
        MONO = 1,
        LOOP = 2,
        NUMBER_OF_CROSS_LINK_TYPES
      };

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

    std::vector<PeptideHit::FragmentAnnotation> frag_annotations;

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


  class OPENMS_DLLAPI OpenProXLUtils
  {

    public:
      struct XLPrecursor
      {
        float precursor_mass;
        unsigned int alpha_index;
        unsigned int beta_index;
      };

      // comparator for sorting XLPrecursor vectors and using upper_bound and lower_bound using only a precursor mass
      struct XLPrecursorComparator
      {
        bool operator() (XLPrecursor& a, XLPrecursor& b) const
        {
          return a.precursor_mass < b.precursor_mass;
        }
        bool operator() (XLPrecursor& a, const double& b) const
        {
          return a.precursor_mass < b;
        }
        bool operator() (const double& a, XLPrecursor& b) const
        {
          return a < b.precursor_mass;
        }
      };


  enum PeptidePosition
  {
    INTERNAL = 0,
    C_TERM = 1,
    N_TERM = 2
  };

  struct AASeqWithMass
  {
    double peptide_mass;
    AASequence peptide_seq;
    PeptidePosition position;
  };

  struct AASeqWithMassComparator
  {
    bool operator() (AASeqWithMass a, AASeqWithMass b) const
    {
      return a.peptide_mass < b.peptide_mass;
    }
    bool operator() (AASeqWithMass a, const double b) const
    {
      return a.peptide_mass < b;
    }
    bool operator() (const double a, AASeqWithMass b) const
    {
      return a < b.peptide_mass;
    }
  };

  struct PreprocessedPairSpectra
  {
    // pre-initialize so we can simply std::swap the spectra (no synchronization in multi-threading context needed as we get no reallocation of the PeakMapreprocessed_pair_spectra.
    PeakMap spectra_common_peaks; // merge spectrum of common peaks (present in both spectra)
    PeakMap spectra_xlink_peaks; // Xlink peaks in the light spectrum (common peaks between spectra_light_different and spectra heavy_to_light)
    PeakMap spectra_all_peaks;

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

    static std::vector<XLPrecursor> enumerateCrossLinksAndMasses_(const std::vector<OpenProXLUtils::AASeqWithMass>&  peptides, double cross_link_mass_light, const DoubleList& cross_link_mass_mono_link, const StringList& cross_link_residue1, const StringList& cross_link_residue2, std::vector< double >& spectrum_precursors, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm);

    static std::vector<ResidueModification> getModificationsFromStringList(StringList modNames);

    static std::vector<OpenProXLUtils::AASeqWithMass> digestDatabase(std::vector<FASTAFile::FASTAEntry> fasta_db, EnzymaticDigestion digestor, Size min_peptide_length, StringList cross_link_residue1, StringList cross_link_residue2, std::vector<ResidueModification> fixed_modifications, std::vector<ResidueModification> variable_modifications, Size max_variable_mods_per_peptide, Size count_proteins = 0, Size count_peptides = 0, bool n_term_linker = false, bool c_term_linker = false);

    static std::vector <ProteinProteinCrossLink> buildCandidates(const std::vector< OpenProXLUtils::XLPrecursor > & candidates, const std::vector<OpenProXLUtils::AASeqWithMass> & peptide_masses, const StringList & cross_link_residue1, const StringList & cross_link_residue2, double cross_link_mass, const DoubleList & cross_link_mass_mono_link, double precursor_mass, double allowed_error, String cross_link_name, bool n_term_linker, bool c_term_linker);

    static void buildFragmentAnnotations(std::vector<PeptideHit::FragmentAnnotation> & frag_annotations, const std::vector< std::pair< Size, Size > > & matching, const PeakSpectrum & theoretical_spectrum, const PeakSpectrum & experiment_spectrum);

    static void buildPeptideIDs(std::vector<PeptideIdentification> & peptide_ids, const std::vector< CrossLinkSpectrumMatch > & top_csms_spectrum, std::vector< std::vector< CrossLinkSpectrumMatch > > & all_top_csms, Size all_top_csms_current_index, const PeakMap & spectra, Size scan_index, Size scan_index_heavy = NULL);

  };

}


#endif










