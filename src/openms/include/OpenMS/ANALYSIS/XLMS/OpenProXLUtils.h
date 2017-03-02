// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
//    PeptideIdentification *peptide_id = NULL;

    bool operator<(const CrossLinkSpectrumMatch& other) const
    {
      return score < other.score;
    }

    bool operator==(const CrossLinkSpectrumMatch& other) const
    {
      return cross_link.alpha == other.cross_link.alpha &&
           cross_link.beta == other.cross_link.beta &&
           cross_link.cross_link_position == other.cross_link.cross_link_position &&
           score == other.score &&
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
           rank == other.rank;
    }

  };


  struct OPENMS_DLLAPI OpenProXLUtils
  {

    struct XLPrecursor
    {
      double precursor_mass;
      Size alpha_index;
      Size beta_index;

      bool operator<(const XLPrecursor& other) const
      {
        return precursor_mass < other.precursor_mass;
      }

      bool operator<(const double other) const
      {
        return precursor_mass < other;
      }
  };

  enum PeptidePosition
  {
    INTERNAL = 0,
    C_TERM = 1,
    N_TERM = 2
  };

  struct PeptideMass
  {
    double peptide_mass;
    AASequence peptide_seq;
    PeptidePosition position;

    bool operator<(const PeptideMass& other) const
    {
      return peptide_mass < other.peptide_mass;
    }

    bool operator<(const double other) const
    {
      return peptide_mass < other;
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

    //bool operator< (const double other, const XLPrecursor& pre) {return other < pre.precursor_mass;}

    static float preScore(Size matchedAlpha, Size ionsAlpha, Size matchedBeta, Size ionsBeta);
    static float preScore(Size matchedAlpha, Size ionsAlpha);

//    static double binomialCoefficient(int n, int k);
    static double cumulativeBinomial(Size n, Size k, double p);

    static std::vector<XLPrecursor> enumerateCrossLinksAndMasses_(const std::vector<OpenProXLUtils::PeptideMass>&  peptides, double cross_link_mass_light, const DoubleList& cross_link_mass_mono_link, const StringList& cross_link_residue1, const StringList& cross_link_residue2, std::vector< double >& spectrum_precursors, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm);


    static double match_odds_score(const PeakSpectrum& theoretical_spec,  const std::vector< std::pair< Size, Size > >& matched_spec, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool is_xlink_spectrum, Size n_charges = 1);

//    template <typename SpectrumType1, typename SpectrumType2>
//    static std::vector< double > xCorrelation(const SpectrumType1 & spec1, const SpectrumType2 & spec2, Int maxshift, double tolerance);

    static double weighted_TIC_score(Size alpha_size, Size beta_size, double intsum_alpha, double intsum_beta, double intsum, double total_current, bool type_is_cross_link);

    // Sum of matched ion intesity, for Intsum score and %TIC score
//    static double matched_current_chain(const std::vector< std::pair< Size, Size > >& matched_spec_common, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks, const PeakSpectrum& spectrum_common_peaks, const PeakSpectrum& spectrum_xlink_peaks);


    static void writeXQuestXML(String out_file, String base_name, const std::vector< PeptideIdentification >& peptide_ids, const std::vector< std::vector< CrossLinkSpectrumMatch > >& all_top_csms, const PeakMap& spectra,
                                                  String precursor_mass_tolerance_unit, String fragment_mass_tolerance_unit, double precursor_mass_tolerance, double fragment_mass_tolerance, double fragment_mass_tolerance_xlinks, String cross_link_name,
                                                  double cross_link_mass_light, DoubleList cross_link_mass_mono_link, String in_fasta, String in_decoy_fasta, StringList cross_link_residue1, StringList cross_link_residue2, double cross_link_mass_iso_shift, String enzyme_name, Size missed_cleavages);

    static void writeXQuestXMLSpec(String out_file, String base_name, const PreprocessedPairSpectra& preprocessed_pair_spectra, const std::vector< std::pair<Size, Size> >& spectrum_pairs, const std::vector< std::vector< CrossLinkSpectrumMatch > >& all_top_csms, const PeakMap& spectra);

    static void writeXQuestXMLSpec(String out_file, String base_name, const std::vector< std::vector< CrossLinkSpectrumMatch > >& all_top_csms, const PeakMap& spectra);

    static PeakSpectrum mergeAnnotatedSpectra(PeakSpectrum & first_spectrum, PeakSpectrum & second_spectrum);

    static void nLargestSpectrumFilter(PeakSpectrum spectrum, int peak_count);

    static void wrap_(const String& input, Size width, String & output);

    static String getxQuestBase64EncodedSpectrum(const PeakSpectrum& spec, String header);

    static std::vector<ResidueModification> getModificationsFromStringList(StringList modNames);

    static void preprocessSpectraLabeled(PeakMap& exp, double fragment_mass_tolerance_xlinks, bool fragment_mass_tolerance_unit_ppm);

    static void getSpectrumAlignment(std::vector<std::pair<Size, Size> > & alignment, const PeakSpectrum & s1, const PeakSpectrum & s2, double tolerance, bool relative_tolerance, double intensity_cutoff = 0.0);

    static PeakSpectrum deisotopeAndSingleChargeMSSpectrum(PeakSpectrum& old_spectrum, Int min_charge, Int max_charge, double fragment_tolerance, bool fragment_tolerance_unit_ppm, bool keep_only_deisotoped = false, Size min_isopeaks = 3, Size max_isopeaks = 10, bool make_single_charged = false);

    static std::vector<OpenProXLUtils::PeptideMass> digestDatabase(std::vector<FASTAFile::FASTAEntry> fasta_db, EnzymaticDigestion digestor, Size min_peptide_length, StringList cross_link_residue1, StringList cross_link_residue2, std::vector<ResidueModification> fixed_modifications, std::vector<ResidueModification> variable_modifications, Size max_variable_mods_per_peptide, Size count_proteins = 0, Size count_peptides = 0, bool n_term_linker = false, bool c_term_linker = false);

    static std::vector <ProteinProteinCrossLink> buildCandidates(const std::vector< OpenProXLUtils::XLPrecursor > & candidates, const std::vector<OpenProXLUtils::PeptideMass> & peptide_masses, const StringList & cross_link_residue1, const StringList & cross_link_residue2, double cross_link_mass, const DoubleList & cross_link_mass_mono_link, double precursor_mass, double allowed_error, String cross_link_name, bool n_term_linker, bool c_term_linker);

    static void buildFragmentAnnotations(std::vector<PeptideHit::FragmentAnnotation> & frag_annotations, const std::vector< std::pair< Size, Size > > & matching, const PeakSpectrum & theoretical_spectrum, const PeakSpectrum & experiment_spectrum);

    static void buildPeptideIDs(std::vector<PeptideIdentification> & peptide_ids, const std::vector< CrossLinkSpectrumMatch > & top_csms_spectrum, std::vector< std::vector< CrossLinkSpectrumMatch > > & all_top_csms, Size all_top_csms_current_index, const PeakMap & spectra, Size scan_index, Size scan_index_heavy = NULL);

    // Sum of matched ion intensity, for Intsum score and %TIC score
    static double matched_current_chain(const std::vector< std::pair< Size, Size > >& matched_spec_common, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks, const PeakSpectrum& spectrum_common_peaks, const PeakSpectrum& spectrum_xlink_peaks);


//    static double total_matched_current(const std::vector< std::pair< Size, Size > >& matched_spec_common_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_common_beta, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_beta, const PeakSpectrum& spectrum_common_peaks, const PeakSpectrum& spectrum_xlink_peaks);

    static double total_matched_current(const std::vector< std::pair< Size, Size > >& matched_spec_common_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_common_beta, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_beta, const PeakSpectrum& spectrum_common_peaks, const PeakSpectrum& spectrum_xlink_peaks);

    // Cross-correlation, with shifting the second spectrum from -maxshift to +maxshift of tolerance bins (Tolerance in Da, a constant binsize)
    static std::vector< double > xCorrelation(const PeakSpectrum & spec1, const PeakSpectrum & spec2, Int maxshift, double tolerance);


  };

  static bool operator< (const double other, const OpenProXLUtils::XLPrecursor& pre) {return other < pre.precursor_mass;}
  static bool operator< (const double other, const OpenProXLUtils::PeptideMass& pre) {return other < pre.peptide_mass;}

}


#endif










