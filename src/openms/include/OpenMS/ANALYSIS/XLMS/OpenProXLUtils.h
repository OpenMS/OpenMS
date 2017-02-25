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
        NUMBER_OF_TERM_SPECIFICITY
      };

      AASequence alpha; // longer peptide
      AASequence beta; // shorter peptide (empty for mono-link), tie bracker: mass then lexicographical
      std::pair<SignedSize, SignedSize> cross_link_position; // index in alpha, beta or between alpha, alpha in loop-links
      double cross_linker_mass;
      String cross_linker_name;
      ResidueModification::Term_Specificity term_spec_alpha;
      ResidueModification::Term_Specificity term_spec_beta;

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




    //bool operator< (const double other, const XLPrecursor& pre) {return other < pre.precursor_mass;}

    static float preScore(Size matchedAlpha, Size ionsAlpha, Size matchedBeta, Size ionsBeta);
    static float preScore(Size matchedAlpha, Size ionsAlpha);

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

    static PeakSpectrum mergeAnnotatedSpectra(PeakSpectrum & first_spectrum, PeakSpectrum & second_spectrum);

    static void nLargestSpectrumFilter(PeakSpectrum spectrum, int peak_count);


    // Sum of matched ion intesity, for Intsum score and %TIC score
    template <typename SpectrumType1>
    static double matched_current_chain(const std::vector< std::pair< Size, Size > >& matched_spec_common, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks, const SpectrumType1& spectrum_common_peaks, const SpectrumType1& spectrum_xlink_peaks)
    {
      double intsum = 0;
      for (SignedSize j = 0; j < static_cast<SignedSize>(matched_spec_common.size()); ++j)
      {
        intsum += spectrum_common_peaks[matched_spec_common[j].second].getIntensity();
      }
      for (SignedSize j = 0; j < static_cast<SignedSize>(matched_spec_xlinks.size()); ++j)
      {
        intsum += spectrum_xlink_peaks[matched_spec_xlinks[j].second].getIntensity();
      }
      return intsum;
    }

//    static double total_matched_current(const std::vector< std::pair< Size, Size > >& matched_spec_common_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_common_beta, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_beta, const PeakSpectrum& spectrum_common_peaks, const PeakSpectrum& spectrum_xlink_peaks);

    template <typename SpectrumType1>
    static double total_matched_current(const std::vector< std::pair< Size, Size > >& matched_spec_common_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_common_beta, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_beta, const SpectrumType1& spectrum_common_peaks, const SpectrumType1& spectrum_xlink_peaks)
    {
      // make vectors of matched peak indices
      double intsum = 0;
      std::vector< Size > indices_common;
      std::vector< Size > indices_xlinks;
      for (Size j = 0; j < matched_spec_common_alpha.size(); ++j)
      {
        indices_common.push_back(matched_spec_common_alpha[j].second);
      }
      for (Size j = 0; j < matched_spec_common_beta.size(); ++j)
      {
        indices_common.push_back(matched_spec_common_beta[j].second);
      }
      for (Size j = 0; j < matched_spec_xlinks_alpha.size(); ++j)
      {
        indices_xlinks.push_back(matched_spec_xlinks_alpha[j].second);
      }
      for (Size j = 0; j < matched_spec_xlinks_beta.size(); ++j)
      {
        indices_xlinks.push_back(matched_spec_xlinks_beta[j].second);
      }

      // make the indices in the vectors unique
      sort(indices_common.begin(), indices_common.end());
      sort(indices_xlinks.begin(), indices_xlinks.end());
      std::vector< Size >::iterator last_unique_common = unique(indices_common.begin(), indices_common.end());
      std::vector< Size >::iterator last_unique_xlinks = unique(indices_xlinks.begin(), indices_xlinks.end());
      indices_common.erase(last_unique_common, indices_common.end());
      indices_xlinks.erase(last_unique_xlinks, indices_xlinks.end());

      // sum over intensities under the unique indices
      for (Size j = 0; j < indices_common.size(); ++j)
      {
        intsum += spectrum_common_peaks[indices_common[j]].getIntensity();
      }
      for (Size j = 0; j < indices_xlinks.size(); ++j)
      {
        intsum += spectrum_xlink_peaks[indices_xlinks[j]].getIntensity();
      }

      return intsum;
    }

    // Cross-correlation, with shifting the second spectrum from -maxshift to +maxshift of tolerance bins (Tolerance in Da, a constant binsize)
    template <typename SpectrumType1, typename SpectrumType2>
    static std::vector< double > xCorrelation(const SpectrumType1 & spec1, const SpectrumType2 & spec2, Int maxshift, double tolerance)
    {
      // generate vector of results, filled with zeroes
      std::vector< double > results(maxshift * 2 + 1, 0);

      // return 0 = no correlation, either positive nor negative, when one of the spectra is empty (e.g. when no common ions or xlink ions could be matched between light and heavy spectra)
      if (spec1.size() == 0 || spec2.size() == 0) {
        return results;
      }

      double maxionsize = std::max(spec1[spec1.size()-1].getMZ(), spec2[spec2.size()-1].getMZ());
      Int table_size = ceil(maxionsize / tolerance)+1;
      std::vector< double > ion_table1(table_size, 0);
      std::vector< double > ion_table2(table_size, 0);

      // Build tables of the same size, each bin has the size of the tolerance
      for (Size i = 0; i < spec1.size(); ++i)
      {
        Size pos = static_cast<Size>(ceil(spec1[i].getMZ() / tolerance));
        // TODO this line leads to using real intensities
  //      ion_table1[pos] = spec1[i].getIntensity();
        // TODO this line leads to using intensities normalized to 10
        ion_table1[pos] = 10.0;
      }
      for (Size i = 0; i < spec2.size(); ++i)
      {
        Size pos =static_cast<Size>(ceil(spec2[i].getMZ() / tolerance));
        // TODO this line leads to using real intensities
  //      ion_table2[pos] = spec2[i].getIntensity();
        // TODO this line leads to using intensities normalized to 10
        ion_table2[pos] = 10.0;
      }

      // Compute means for real intensities
      double mean1 = (std::accumulate(ion_table1.begin(), ion_table1.end(), 0.0)) / table_size;
      double mean2 = (std::accumulate(ion_table2.begin(), ion_table2.end(), 0.0)) / table_size;

      // Compute denominator
      double s1 = 0;
      double s2 = 0;
      for (Int i = 0; i < table_size; ++i)
      {
        s1 += pow((ion_table1[i] - mean1), 2);
        s2 += pow((ion_table2[i] - mean2), 2);
      }
      double denom = sqrt(s1 * s2);

      // Calculate correlation for each shift
      for (Int shift = -maxshift; shift <= maxshift; ++shift)
      {
        double s = 0;
        for (Int i = 0; i < table_size; ++i)
        {
          Int j = i + shift;
          if ( (j >= 0) && (j < table_size))
          {
            s += (ion_table1[i] - mean1) * (ion_table2[j] - mean2);
          }
        }
        if (denom > 0)
        {
          results[shift + maxshift] = s / denom;
        }
      }
      return results;
    }

  };

  static bool operator< (const double other, const OpenProXLUtils::XLPrecursor& pre) {return other < pre.precursor_mass;}
  static bool operator< (const double other, const OpenProXLUtils::PeptideMass& pre) {return other < pre.peptide_mass;}

}


#endif










