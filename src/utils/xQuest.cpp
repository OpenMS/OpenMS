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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/EnzymesDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>

// preprocessing and filtering
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLinks.h>

#include <iostream>
#include <cmath>

using namespace std;
using namespace OpenMS;


#ifdef _OPENMP
#include <omp.h>
#define NUMBER_OF_THREADS (omp_get_num_threads())
#else
#define NUMBER_OF_THREADS (1)
#endif


/**
    @page UTILS_xQuest xQuest

    @brief Perform protein-protein cross-linking experiments.

    <CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ xQuest \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
        </tr>
    </table>
</CENTER>

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_xQuest.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_xQuest.html
*/

class TOPPxQuest :
  public TOPPBase
{
public:
  TOPPxQuest() :
    TOPPBase("xQuest", "Tool for protein-protein cross linking using the xQuest algorithm.", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    // input files
    registerInputFile_("in", "<file>", "", "Input file containing the spectra.");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("consensus", "<file>", "", "Input file containing the linked mass peaks.");
    setValidFormats_("consensus", ListUtils::create<String>("consensusXML"));

    registerInputFile_("database", "<file>", "", "Input file containing the protein database.");
    setValidFormats_("database", ListUtils::create<String>("fasta"));

    registerTOPPSubsection_("precursor", "Precursor (Parent Ion) Options");
    registerDoubleOption_("precursor:mass_tolerance", "<tolerance>", 10.0, "Width of precursor mass tolerance window", false);

    StringList precursor_mass_tolerance_unit_valid_strings;
    precursor_mass_tolerance_unit_valid_strings.push_back("ppm");
    precursor_mass_tolerance_unit_valid_strings.push_back("Da");

    registerStringOption_("precursor:mass_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance.", false, false);
    setValidStrings_("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);

    registerIntOption_("precursor:min_charge", "<num>", 3, "Minimum precursor charge to be considered.", false, true);
    registerIntOption_("precursor:max_charge", "<num>", 6, "Maximum precursor charge to be considered.", false, true);

    registerTOPPSubsection_("fragment", "Fragments (Product Ion) Options");
    registerDoubleOption_("fragment:mass_tolerance", "<tolerance>", 0.2, "Fragment mass tolerance", false);

    StringList fragment_mass_tolerance_unit_valid_strings;
    fragment_mass_tolerance_unit_valid_strings.push_back("ppm");
    fragment_mass_tolerance_unit_valid_strings.push_back("Da");

    registerStringOption_("fragment:mass_tolerance_unit", "<unit>", "Da", "Unit of fragment m", false, false);
    setValidStrings_("fragment:mass_tolerance_unit", fragment_mass_tolerance_unit_valid_strings);

    registerTOPPSubsection_("modifications", "Modifications Options");
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("modifications:fixed", "<mods>", ListUtils::create<String>(""), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'", false);
    setValidStrings_("modifications:fixed", all_mods);
    registerStringList_("modifications:variable", "<mods>", ListUtils::create<String>(""), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'", false);
    setValidStrings_("modifications:variable", all_mods);
    registerIntOption_("modifications:variable_max_per_peptide", "<num>", 2, "Maximum number of residues carrying a variable modification per candidate peptide", false, false);

    registerTOPPSubsection_("peptide", "Peptide Options");
    registerIntOption_("peptide:min_size", "<num>", 7, "Minimum size a peptide must have after digestion to be considered in the search.", false, true);
    registerIntOption_("peptide:missed_cleavages", "<num>", 1, "Number of missed cleavages.", false, false);
    vector<String> all_enzymes;
    EnzymesDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("peptide:enzyme", "<cleavage site>", "Trypsin", "The enzyme used for peptide digestion.", false);
    setValidStrings_("peptide:enzyme", all_enzymes);


    registerTOPPSubsection_("cross_linker", "Cross Linker Options");
    registerDoubleOption_("cross_linker:mass_light", "<mass>", 156.078644, "Mass of the light cross-linker", false);
    registerDoubleOption_("cross_linker:mass_heavy", "<mass>", 168.153965, "Mass of the heavy cross-linker", false);
    registerDoubleOption_("cross_linker:mass_loss_type2", "<mass>", 18.01056, "Mass difference observed in an intra or inter peptide link", false);

    registerTOPPSubsection_("algorithm", "Algorithm Options");
    registerStringOption_("algorithm:candidate_search", "<param>", "index", "Mode used to generate candidate peptides.", false, false);
    StringList candidate_search_modes_strings;
    candidate_search_modes_strings.push_back("index");
    candidate_search_modes_strings.push_back("enumeration");
    setValidStrings_("algorithm:candidate_search", candidate_search_modes_strings);

    // output file
    registerOutputFile_("out", "<file>", "", "Result file\n");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
  }

  vector<ResidueModification> getModifications_(StringList modNames)
  {
    vector<ResidueModification> modifications;

    // iterate over modification names and add to vector
    for (StringList::iterator mod_it = modNames.begin(); mod_it != modNames.end(); ++mod_it)
    {
      String modification(*mod_it);
      modifications.push_back(ModificationsDB::getInstance()->getModification(modification));
    }

    return modifications;
  }

  // check if for minimum size
  class HasInvalidPeptideLengthPredicate
  {
      public:
        explicit HasInvalidPeptideLengthPredicate(Size min_size)
          :min_size_(min_size)
        {
        }

        bool operator()(const AASequence& aas)
        {
          return (aas.size() < min_size_);
        }
    private:
        Size min_size_;
  };



  void preprocessSpectra_(PeakMap& exp)
  {
    // filter MS2 map
    // remove 0 intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakMap(exp);
 
    Normalizer normalizer;
    normalizer.filterPeakMap(exp);

    // sort by rt
    exp.sortSpectra(false);

    // filter settings
    WindowMower window_mower_filter;
    Param filter_param = window_mower_filter.getParameters();
    filter_param.setValue("windowsize", 100.0, "The size of the sliding window along the m/z axis.");
    filter_param.setValue("peakcount", 20, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    window_mower_filter.setParameters(filter_param);
    NLargest nlargest_filter = NLargest(400);
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size(); ++exp_index)
    {
      // sort by mz
      exp[exp_index].sortByPosition();
  
      // remove noise
      window_mower_filter.filterPeakSpectrum(exp[exp_index]);
      nlargest_filter.filterPeakSpectrum(exp[exp_index]);
  
      // sort (nlargest changes order)
      exp[exp_index].sortByPosition();
     }
   }

  // xQuest, fast pre-Score for x-links (type 2)
  // required: numbers of peaks for each chain, and how many of them were matched
  float preScore(Size matchedAlpha, Size ionsAlpha, Size matchedBeta, Size ionsBeta)
  {
    if ( (ionsAlpha > 0) && (ionsBeta > 0) )
    {
      float result = sqrt(((float) matchedAlpha / (float) ionsAlpha) * ((float) matchedBeta / (float) ionsBeta));
      return result;
    } else
    {
      return 0.0;
    }
  }

  // xQuest, fast pre-Score for Mono links and Loop links (type 0 and type 1)
  float preScore(Size matchedAlpha, Size ionsAlpha)
  {
    if (ionsAlpha > 0)
    {
      float result = (float) matchedAlpha / (float) ionsAlpha;
      return result;
    } else
    {
      return 0.0;
    }
  }

  // Statistics/Combinatorics functions for match-odds score
  // Standard choose n over k
  Int choose(Size n, Size k)
  {
        Size result = 1;
        Size j = 1;
        if ( k > n || k < 0) return 0;

	while (j <= k)
	{
		result *= n--;
		result /= j++;
	}
	return result;
  }

  // Standard cumulative binomial distribution
  double cumulativeBinomial(Size n, Size k, double p)
  {
        double p_cumul = 0.0;
        if (p < 1e-99) return static_cast<double>(k == 0); //  (not true/false, but 1/0 as probability)
        if (1 - p < 1e-99) return static_cast<double>(k != n); //
        if (k > n)  return 1.0;

	for (Size j = 0; j < k; j++)
	{
		p_cumul += choose(n, j) * pow(p,  j) * pow((1-p), (n-j));
	}
	return p_cumul;
  }


  // A priori probability of a random match given info about the theoretical spectrum, (range = maxmz - minmz)
  double aPrioriProb(double ms2_tolerance, Size n_peaks, double range, Size n_charges)
  {
    double a_priori_probability= (1 - ( pow( (1 - 2 * ms2_tolerance/ (0.5 * range)),  (n_peaks / n_charges))));
    return a_priori_probability;
  }
     

  struct PreprocessedPairSpectra_
  {
    // pre-initialize so we can simply std::swap the spectra (no synchronization in multi-threading context needed as we get no reallocation of the PeakMaps) 
    PeakMap spectra_light_different; // peaks in light spectrum after common peaks have been removed
    PeakMap spectra_heavy_different; // peaks in heavy spectrum after common peaks have been removed
    PeakMap spectra_heavy_to_light; // heavy peaks transformed to light ones and after common peaks have been removed
    PeakMap spectra_common_peaks; // merge spectrum of common peaks (present in both spectra)
    PeakMap spectra_xlink_peaks; // Xlink peaks in the light spectrum (common peaks between spectra_light_different and spectra heavy_to_light)

    PreprocessedPairSpectra_(Size size)
    {
      for (Size i = 0; i != size; ++i)
      {
        spectra_light_different.addSpectrum(PeakSpectrum());
        spectra_heavy_different.addSpectrum(PeakSpectrum());
        spectra_heavy_to_light.addSpectrum(PeakSpectrum());
        spectra_common_peaks.addSpectrum(PeakSpectrum());
        spectra_xlink_peaks.addSpectrum(PeakSpectrum());
      }
    }  
  };


  // create common / shifted peak spectra for all pairs
  PreprocessedPairSpectra_ preprocessPairs_(const PeakMap& spectra, const map<Size, Size>& map_light_to_heavy, const SpectrumAlignment& ms2_aligner, const double cross_link_mass_light, const double cross_link_mass_heavy)
  {
    PreprocessedPairSpectra_ ps(spectra.size());
 
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize scan_index = 0; scan_index < (SignedSize)spectra.size(); ++scan_index)
    {
      // assume that current MS2 corresponds to a light spectrum
      const PeakSpectrum& spectrum_light = spectra[scan_index];

      // map light to heavy
      map<Size, Size>::const_iterator scan_index_light_it = map_light_to_heavy.find(scan_index);

      if (scan_index_light_it != map_light_to_heavy.end())
      {
        const Size scan_index_heavy = scan_index_light_it->second;
        const PeakSpectrum& spectrum_heavy = spectra[scan_index_heavy];
//      cout << "light spectrum index: " << scan_index << " heavy spectrum index: " << scan_index_heavy << endl;
        std::vector< std::pair< Size, Size > > matched_fragments_without_shift;
        ms2_aligner.getSpectrumAlignment(matched_fragments_without_shift, spectrum_light, spectrum_heavy);

        // different fragments may carry light or heavy cross-linker.             
        PeakSpectrum spectrum_heavy_different;
        PeakSpectrum spectrum_light_different;

        // TODO: maybe speed this up - can be done in linear time
        for (Size i = 0; i != spectrum_light.size(); ++i)
        {
          bool found = false;
          for (Size j = 0; j != matched_fragments_without_shift.size(); ++j)
          {
            if (matched_fragments_without_shift[j].first == i) { found = true; break; }
          }
          if (!found)
          {
            spectrum_light_different.push_back(spectrum_light[i]);
          }
        }
        for (Size i = 0; i != spectrum_heavy.size(); ++i)
        {
          bool found = false;
          for (Size j = 0; j != matched_fragments_without_shift.size(); ++j)
          {
            if (matched_fragments_without_shift[j].second == i) { found = true; break; }
          }
          if (!found)
          {
            spectrum_heavy_different.push_back(spectrum_heavy[i]);
          }
        }
        std::swap(ps.spectra_light_different[scan_index], spectrum_light_different);
        std::swap(ps.spectra_heavy_different[scan_index], spectrum_heavy_different);

        // transform by m/z difference between unlabeled and labeled cross-link to make heavy and light comparable.
        // TODO: important: assume different charged MS2 fragments. Now only single charged ones are assumed

        // transform heavy spectrum to light spectrum
        PeakSpectrum spectrum_heavy_to_light;
        for (Size i = 0; i != spectrum_heavy_different.size(); ++i)
        {
          Peak1D p = spectrum_heavy_different[i];
          p.setMZ(p.getMZ() - (cross_link_mass_heavy - cross_link_mass_light));
          spectrum_heavy_to_light.push_back(p); 
        }
        std::swap(ps.spectra_heavy_to_light[scan_index], spectrum_heavy_to_light);

        /*
        // transform heavy spectrum to artifical spectrum with cross-linker completely removed
        PeakSpectrum spectrum_heavy_no_linker;
        for (Size i = 0; i != spectrum_heavy_different.size(); ++i)
        {
          Peak1D p = spectrum_heavy_different[i];
          p.setMZ(p.getMZ() - cross_link_mass_heavy);
          spectrum_heavy_no_linker.push_back(p); 
        }

        // transform light spectrum to artifical spectrum with cross-linker completely removed
        PeakSpectrum spectrum_light_no_linker;
        for (Size i = 0; i != spectrum_light_different.size(); ++i)
        {
          Peak1D p = spectrum_light_different[i];
          p.setMZ(p.getMZ() - cross_link_mass_light);
          spectrum_light_no_linker.push_back(p); 
        }
        */

        // align potentially shifted peaks from light MS2 with potentially shifted peaks from heavy (after transformation to resemble the light MS2)
        // matching fragments are potentially carrying the cross-linker
        std::vector< std::pair< Size, Size > > matched_fragments_with_shift;
        ms2_aligner.getSpectrumAlignment(matched_fragments_with_shift, spectrum_light_different, spectrum_heavy_to_light);
        PeakSpectrum xlink_peaks;
        for (Size i = 0; i != matched_fragments_with_shift.size(); ++i)
        {
          xlink_peaks.push_back(spectrum_light[matched_fragments_with_shift[i].first]);
        }
        xlink_peaks.sortByPosition();
        std::swap(ps.spectra_xlink_peaks[scan_index], xlink_peaks);

#ifdef DEBUG_XQUEST
        cout << "Common peaks: " << matched_fragments_without_shift.size() << " different peaks: " << spectrum_light.size() - matched_fragments_without_shift.size() << ", " << spectrum_heavy.size() - matched_fragments_without_shift.size() << endl;
        cout << "Matched shifted peaks: " << matched_fragments_with_shift.size() << " unexplained peaks: " << spectrum_light_different.size() - matched_fragments_with_shift.size() << ", " << spectrum_heavy_to_light.size() - matched_fragments_with_shift.size() << endl;
#endif
        // generate common peaks spectrum TODO: check if only light / only heavy or e.g. mean of both peaks are used here
        PeakSpectrum common_peaks;
        for (Size i = 0; i != matched_fragments_without_shift.size(); ++i)
        {
          common_peaks.push_back(spectrum_light[matched_fragments_without_shift[i].first]);
        }
        common_peaks.sortByPosition();
#ifdef DEBUG_XQUEST
        cout << "Peaks to match: " << common_peaks.size() << endl;
#endif
        std::swap(ps.spectra_common_peaks[scan_index], common_peaks);
      }
    }
    return ps;
  }

  struct OPENMS_DLLAPI MatchedIonCount
  {
    // Note: code is optimize for theo_spectrum to contain less peaks than the exp_spectrum
    static Size compute(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const MSSpectrum<RichPeak1D>& theo_spectrum, const MSSpectrum<Peak1D>& exp_spectrum)
    {
      Size matches(0), i(0), j(0);
      const Size count_exp_spectrum = exp_spectrum.size();
      const Size count_theo_spectrum = theo_spectrum.size();
      do
      { 
        // advance i until in or right of tolerance window
        while (i < count_exp_spectrum && j < count_theo_spectrum)
        {
          const double tolerance_Th = fragment_mass_tolerance_unit_ppm ? theo_spectrum[j].getMZ() * 1e-6 * fragment_mass_tolerance : fragment_mass_tolerance;
          if (exp_spectrum[i].getMZ() < (theo_spectrum[j].getMZ() - tolerance_Th)) 
          {
            ++i;
          }
          else
          {
            break;
          }
        }

        if (i == count_exp_spectrum || j == count_theo_spectrum) return matches;

        // check if in tolerance window (could also be to the right of it) and count as match
        const double tolerance_Th = fragment_mass_tolerance_unit_ppm ? theo_spectrum[j].getMZ() * 1e-6 * fragment_mass_tolerance : fragment_mass_tolerance;
        if (exp_spectrum[i].getMZ() < (theo_spectrum[j].getMZ() + tolerance_Th))
        {
          // matched i and j (can't match again)
          ++matches; 
          ++j; 
          ++i; 
        } 

        if (i == count_exp_spectrum || j == count_theo_spectrum) return matches;

        ++j;
      } while (true);  
      return matches;
    }
  };

  struct OPENMS_DLLAPI HashGrid1D
  {
    // bucket size should be 2.0 * fragment mass tolerance
    // this ensures that two neighboring buckets cover all position possible given the fragment mass tolerance
    HashGrid1D(double min, double max, double bucket_size) : 
      min_(min), 
      max_(max), 
      bucket_size_(bucket_size)
    {
      Size n_buckets = ceil((max - min) / bucket_size) + 1;
      h_.resize(n_buckets);
    }

    void insert(double position, AASequence* v)
    {
      if (position < min_) position = min_;
      if (position > max_) position = max_;

      double bucket_index = (position - min_) / bucket_size_;

      if (bucket_index - (Size)bucket_index <= 0.5)
      {
        h_[bucket_index].push_back(v);
        if (bucket_index >= 0) h_[bucket_index - 1].push_back(v);
      } 
      else
      {
        h_[bucket_index].push_back(v);
        if (bucket_index < h_.size() - 1) h_[bucket_index + 1].push_back(v);
      }
    }

    vector<AASequence*>& get(double position) 
    {
      double bucket_index = (position - min_) / bucket_size_;
      return h_[bucket_index];
    }

    std::vector<vector<AASequence*> > h_; 
    double min_;
    double max_;
    double bucket_size_;
  };


  static multimap<double, pair<AASequence, AASequence> > enumerateCrossLinksAndMasses_(const multimap<StringView, AASequence>&  peptides, double cross_link_mass_light, double cross_link_mass_loss_type2)
  {
    multimap<double, pair<AASequence, AASequence> > mass_to_candidates;
    for (map<StringView, AASequence>::const_iterator a = peptides.begin(); a != peptides.end(); ++a)
    {
      if (a->second.toString().find("K") >= a->second.size()-1)
      {
        continue;
      }

      for (map<StringView, AASequence>::const_iterator b = a; b != peptides.end(); ++b)
      {

        if (b->second.toString().find("K") >= b->second.size()-1)
        {
          continue;
        }
/*
        // Find all positions of lysine (K) in the peptides (possible scross-linking sites)
        vector <Size> K_pos_a;
        vector <Size> K_pos_b;
        for (Size i = 0; i < a->second.size()-1; ++i)
        {
          if (a->second.toString()[i] == 'K') K_pos_a.push_back(i); // TODO: query by StringView (might need operator[] in StringView)
        }
        for (Size i = 0; i < b->second.size()-1; ++i)
        {
          if (b->second.toString()[i] == 'K') K_pos_b.push_back(i);
        }
*/
        // mass peptide1 + mass peptide2 + cross linker mass - cross link loss
        double cross_link_mass = a->second.getMonoWeight() + b->second.getMonoWeight() + cross_link_mass_light - cross_link_mass_loss_type2;
        mass_to_candidates.insert(make_pair(cross_link_mass, make_pair<AASequence, AASequence>(a->second, b->second)));
      }
    }

    return mass_to_candidates;
  }

  ExitCodes main_(int, const char**)
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);

    const string in_mzml(getStringOption_("in"));
    const string in_fasta(getStringOption_("database"));
    const string in_consensus(getStringOption_("consensus"));
    const string out_idxml(getStringOption_("out"));

    Int min_precursor_charge = getIntOption_("precursor:min_charge");
    Int max_precursor_charge = getIntOption_("precursor:max_charge");
    double precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
    bool precursor_mass_tolerance_unit_ppm = (getStringOption_("precursor:mass_tolerance_unit") == "ppm");

    double fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
    bool fragment_mass_tolerance_unit_ppm = (getStringOption_("fragment:mass_tolerance_unit") == "ppm");

    SpectrumAlignment ms2_aligner;
    Param ms2_alinger_param = ms2_aligner.getParameters();
    String ms2_relative_tolerance = fragment_mass_tolerance_unit_ppm ? "true" : "false";
    ms2_alinger_param.setValue("is_relative_tolerance", ms2_relative_tolerance);
    ms2_alinger_param.setValue("tolerance", fragment_mass_tolerance);
    ms2_aligner.setParameters(ms2_alinger_param);

    double cross_link_mass_light = getDoubleOption_("cross_linker:mass_light");
    double cross_link_mass_heavy = getDoubleOption_("cross_linker:mass_heavy");
    double cross_link_mass_loss_type2 = getDoubleOption_("cross_linker:mass_loss_type2");

    StringList fixedModNames = getStringList_("modifications:fixed");
    set<String> fixed_unique(fixedModNames.begin(), fixedModNames.end());

    Size peptide_min_size = getIntOption_("peptide:min_size");

    bool ion_index_mode = (getStringOption_("algorithm:candidate_search") == "index");

    if (fixed_unique.size() != fixedModNames.size())
    {
      cout << "duplicate fixed modification provided." << endl;
      return ILLEGAL_PARAMETERS;
    }

    StringList varModNames = getStringList_("modifications:variable");
    set<String> var_unique(varModNames.begin(), varModNames.end());
    if (var_unique.size() != varModNames.size())
    {
      cout << "duplicate variable modification provided." << endl;
      return ILLEGAL_PARAMETERS;
    }

    vector<ResidueModification> fixed_modifications = getModifications_(fixedModNames);
    vector<ResidueModification> variable_modifications = getModifications_(varModNames);
    Size max_variable_mods_per_peptide = getIntOption_("modifications:variable_max_per_peptide");
    
    // load MS2 map
    PeakMap spectra;
    MzMLFile f;
    f.setLogType(log_type_);

    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(2);
    f.getOptions() = options;
    f.load(in_mzml, spectra);
    spectra.sortSpectra(true);

    // filter noise
    progresslogger.startProgress(0, 1, "Filtering spectra...");
    preprocessSpectra_(spectra);
    progresslogger.endProgress();

    // load linked features
    ConsensusMap cfeatures;
    ConsensusXMLFile cf;
    cf.load(in_consensus, cfeatures); 

    // load fasta database    
    progresslogger.startProgress(0, 1, "Load database from FASTA file...");
    FASTAFile fastaFile;
    vector<FASTAFile::FASTAEntry> fasta_db;
    fastaFile.load(in_fasta, fasta_db);
    progresslogger.endProgress();
    
    const Size missed_cleavages = getIntOption_("peptide:missed_cleavages");
    EnzymaticDigestion digestor;
    digestor.setEnzyme(getStringOption_("peptide:enzyme"));
    digestor.setMissedCleavages(missed_cleavages);

    progresslogger.startProgress(0, (Size)(fasta_db.end() - fasta_db.begin()), "Scoring peptide models against spectra...");

    // lookup for processed peptides. must be defined outside of omp section and synchronized
    multimap<StringView, AASequence> processed_peptides;
    
    // set minimum size of peptide after digestion
    Size min_peptide_length = getIntOption_("peptide:min_size");

    // build multimap of precursor mass to scan index
    multimap<double, Size> multimap_mass_2_scan_index;

    vector<PeptideIdentification> pseudo_ids; // used to map precursor positions to consensus features
    for (PeakMap::ConstIterator s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
    {
      int scan_index = s_it - spectra.begin();
      vector<Precursor> precursor = s_it->getPrecursors();
    
      // there should only one precursor and MS2 should contain at least a few peaks to be considered (e.g. at least for every AA in the peptide)
      if (precursor.size() == 1 && s_it->size() >= peptide_min_size)
      {
        int precursor_charge = precursor[0].getCharge();
        if (precursor_charge < min_precursor_charge || precursor_charge > max_precursor_charge)
        {
          continue;
        }
    
        double precursor_mz = precursor[0].getMZ();
        double precursor_mass = (double) precursor_charge * precursor_mz - (double) precursor_charge * Constants::PROTON_MASS_U;
    
        multimap_mass_2_scan_index.insert(make_pair(precursor_mass, scan_index));
        PeptideIdentification temp_pi;
        temp_pi.setRT(s_it->getRT());
        temp_pi.setMZ(precursor_mz);
        temp_pi.setMetaValue("scan_index", scan_index);
        vector<PeptideHit> temp_hits;
        PeptideHit temp_ph;
        temp_ph.setCharge(precursor_charge);
        temp_hits.push_back(temp_ph);
        temp_pi.setHits(temp_hits);
        pseudo_ids.push_back(temp_pi);
      }
    }

    IDMapper idmapper;
    Param p = idmapper.getParameters();
    p.setValue("rt_tolerance", 30.0);
    p.setValue("mz_tolerance", precursor_mass_tolerance);
    String mz_measure = precursor_mass_tolerance_unit_ppm ? "ppm" : "Da";
    p.setValue("mz_measure", mz_measure);
    p.setValue("mz_reference", "precursor");
    p.setValue("ignore_charge", "false");
    idmapper.setParameters(p);

    vector<ProteinIdentification> protein_ids;
    // protein identifications (leave as is...)
    protein_ids = vector<ProteinIdentification>(1);
    protein_ids[0].setDateTime(DateTime::now());
    protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());
    idmapper.annotate(cfeatures, pseudo_ids, protein_ids, true, true);

    // maps the index of a light precursor peptide to its corresponding heavier partner
    map<Size, Size> map_light_to_heavy;
    for (ConsensusMap::const_iterator cit = cfeatures.begin(); cit != cfeatures.end(); ++cit)
    {
      if (cit->getFeatures().size() == 2 && cit->getPeptideIdentifications().size() == 2)
      {
        const PeptideIdentification& pi_0 = cit->getPeptideIdentifications()[0];
        const PeptideIdentification& pi_1 = cit->getPeptideIdentifications()[1];
        map_light_to_heavy[pi_0.getMetaValue("scan_index")] = pi_1.getMetaValue("scan_index");
      }
    }
   
    cout << "Number of MS2 pairs connceted by consensus feature: " << map_light_to_heavy.size() << endl;

    // create common peak / shifted peak spectra for all pairs
    PreprocessedPairSpectra_ ps = preprocessPairs_(spectra, map_light_to_heavy, ms2_aligner, cross_link_mass_light, cross_link_mass_heavy);
 
    Size count_proteins = 0;
    Size count_peptides = 0;
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize fasta_index = 0; fasta_index < (SignedSize)fasta_db.size(); ++fasta_index)
    {
#ifdef _OPENMP
#pragma omp atomic
#endif
      ++count_proteins;

      IF_MASTERTHREAD
      {
        progresslogger.setProgress((SignedSize)fasta_index * NUMBER_OF_THREADS);
      }

      // store vector of substrings pointing in fasta database (bounded by pairs of begin, end iterators)    
      vector<StringView> current_digest;
      digestor.digestUnmodifiedString(fasta_db[fasta_index].sequence, current_digest, min_peptide_length);

      for (vector<StringView>::iterator cit = current_digest.begin(); cit != current_digest.end(); ++cit)
      {
        bool already_processed = false;
#ifdef _OPENMP
#pragma omp critical (processed_peptides_access)
#endif
        {
          if (processed_peptides.find(*cit) != processed_peptides.end())
          {
            // peptide (and all modified variants) already processed so skip it
            already_processed = true;
          }
        }

        if (already_processed)
        {
          continue;
        }
        if (cit->getString().find('K') >= cit->getString().size()-1)
        {
          continue;
        }
            
#ifdef _OPENMP
#pragma omp atomic
#endif
        ++count_peptides;

        vector<AASequence> all_modified_peptides;

        // generate all modified variants of a peptide
        // Note: no critial section is needed despite ResidueDB not beeing thread sage.
        //       It is only written to on introduction of novel modified residues. These residues have been already added above (single thread context).
        {
          AASequence aas = AASequence::fromString(cit->getString());
          ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications.begin(), fixed_modifications.end(), aas);
          ModifiedPeptideGenerator::applyVariableModifications(variable_modifications.begin(), variable_modifications.end(), aas, max_variable_mods_per_peptide, all_modified_peptides);
        }
        
        for (SignedSize mod_pep_idx = 0; mod_pep_idx < (SignedSize)all_modified_peptides.size(); ++mod_pep_idx)
        {
          const AASequence& candidate = all_modified_peptides[mod_pep_idx];

#ifdef _OPENMP
#pragma omp critical (processed_peptides_access)
#endif
          {
            processed_peptides.insert(pair<StringView, AASequence>(*cit, candidate));
          }
        }
      }
    }

    cout << "Peptide " << processed_peptides.size() << " candidates." << endl;

    //TODO refactor, so that only the used mode is initialized and the pre-scoring code only appears once
    // Initialize enumeration mode
    multimap<double, pair<AASequence, AASequence> > enumerated_cross_link_masses;

     // Initialize Ion Index mode
     //Size counter(0);
     //TODO remove, adapt to ppm
     HashGrid1D hg(0.0, 5000.0, fragment_mass_tolerance);

    if (!ion_index_mode)
    {
      enumerated_cross_link_masses = enumerateCrossLinksAndMasses_(processed_peptides, cross_link_mass_light, cross_link_mass_loss_type2);
            cout << "Enumerated cross-links: " << enumerated_cross_link_masses.size() << endl;
    } else {


        // create spectrum generator
        TheoreticalSpectrumGeneratorXLinks spectrum_generator;

        Size has_aligned_peaks(0);
        Size no_aligned_peaks(0);

        // filtering peptide candidates
        for (map<StringView, AASequence>::iterator a = processed_peptides.begin(); a != processed_peptides.end(); ++a)
        {
          //create theoretical spectrum
          MSSpectrum<RichPeak1D> theo_spectrum = MSSpectrum<RichPeak1D>();

          const AASequence& seq = a->second;
          //add peaks for b and y ions with charge 1
          //cout << a->first.getString() << ":" << a->second.toString() << endl;
          spectrum_generator.getCommonIonSpectrum(theo_spectrum, seq, 1);

          //sort by mz
          theo_spectrum.sortByPosition();

          for (Size i = 0; i != theo_spectrum.size(); ++i)
          {
            hg.insert(theo_spectrum[i].getMZ(), &(a->second)); // TODO add real index here
          }
        }
      }


    // TODO test variable, can be removed
    float pScoreMax =0;
    double TICMax = 0;
    double matchOddsMax = 0;

    // iterate over all spectra
    for (SignedSize scan_index = 0; scan_index < (SignedSize)spectra.size(); ++scan_index)
    {
      // assume that current MS2 corresponds to a light spectrum
      const PeakSpectrum& spectrum_light = spectra[scan_index];
      const double precursor_charge = spectrum_light.getPrecursors()[0].getCharge();
      const double precursor_mz = spectrum_light.getPrecursors()[0].getMZ();
      const double precursor_mass = precursor_mz * static_cast<double>(precursor_charge) - static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U;

      // map light to heavy
      map<Size, Size>::const_iterator scan_index_light_it = map_light_to_heavy.find(scan_index);
   
      if (scan_index_light_it != map_light_to_heavy.end())
      { 
        const Size scan_index_heavy = scan_index_light_it->second;
        const PeakSpectrum& spectrum_heavy = spectra[scan_index_heavy];

        cout << "Pair: " << scan_index << ", " << scan_index_heavy << " (mz, charge, mass) " << precursor_mz << "," << precursor_charge << "," << precursor_mass <<  endl;
//          cout << "light spectrum index: " << scan_index_light << " heavy spectrum index: " << scan_index_heavy << endl;

        // Matching of common peaks (done with 0.2 Da tolerance in xQuest)
//        std::vector< std::pair< Size, Size > > common_peaks; // TODO: use precalculated
//        ms2_aligner.getSpectrumAlignment(common_peaks, spectrum_light, spectrum_heavy);
        PeakSpectrum& common_peaks = ps.spectra_common_peaks[scan_index];

/*
        cout << "Pair: " << a->second << "(" << a->second.getMonoWeight() << ")" << ", " 
             << b->second << "(" << b->second.getMonoWeight() << ") matched to light spectrum " << scan_index_light << " with m/z: " << spectrum_light.getPrecursors()[0].getMZ() << " cross_link_mass: " <<  cross_link_mass << endl;
        cout << "Common peaks: " << common_peaks.size() << " remaining preaks: " << spectrum_light.size() - common_peaks.size() << ", " << spectrum_heavy.size() - common_peaks.size() << endl;
       
*/

        if(common_peaks.size() > 3)
        {
          // determine candidates
          if (ion_index_mode)
          {
            // TODO Use 50 most intense common peaks of exp. spectrum, consider all peptides that produce any of these as theor. common ions
            NLargest nlargest_filter = NLargest(50);
            double most_intensive_peak_mz(0);
            double most_intensive_peak_int(-1);
            for (Size i = 0; i != common_peaks.size(); ++i)
            {
              double current_intensity = common_peaks[i].getIntensity();
              double current_mz = common_peaks[i].getMZ();
              if (current_intensity > most_intensive_peak_int)
              {
                most_intensive_peak_int = current_intensity;
                most_intensive_peak_mz = current_mz;
              }
            }
            const vector<AASequence*> ion_tag_candidates = hg.get(most_intensive_peak_mz);
            cout << "Ion tag candidates before mass filtering: " << ion_tag_candidates.size() << endl;

            vector<pair<AASequence, AASequence> > filtered_candidates;
            for (Size i = 0; i != ion_tag_candidates.size(); ++i)
            {
              const AASequence* peptide_a = ion_tag_candidates[i];
              
              if (peptide_a->toString().find("K") >= peptide_a->size()-1)
              {
                continue;
              }

              for (Size j = i + 1; j < ion_tag_candidates.size(); ++j)
              {
                const AASequence* peptide_b = ion_tag_candidates[j];

                if (peptide_b->toString().find("K") >= peptide_b->size()-1)
                {
                  continue;
                }

                double cross_link_mass = peptide_a->getMonoWeight() + peptide_b->getMonoWeight() + cross_link_mass_light - cross_link_mass_loss_type2;
                double error_Da = abs(cross_link_mass - precursor_mass);
                if (error_Da < precursor_mass_tolerance)
                {
                  filtered_candidates.push_back(make_pair(*peptide_a, *peptide_b));
                }                
              } 
            } 
            cout << "Ion tag candidates after mass filtering: " << filtered_candidates.size() << endl;

            for (Size i = 0; i != filtered_candidates.size(); ++i) 
            {
              const pair<AASequence, AASequence>& candidate = filtered_candidates[i];
	      TheoreticalSpectrumGeneratorXLinks specGen;
              RichPeakSpectrum theoretical_spec;
              TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink cross_link_candidate;
              cross_link_candidate.alpha = candidate.first;
              cross_link_candidate.beta = candidate.second;
              cross_link_candidate.cross_link_position.first = candidate.first.toString().find('K');
              cross_link_candidate.cross_link_position.second = candidate.second.toString().find('K');

              specGen.getSpectrum(theoretical_spec, cross_link_candidate, 1);

              std::vector< std::pair< Size, Size > > matched_fragments_theor_spec;
              ms2_aligner.getSpectrumAlignment(matched_fragments_theor_spec, theoretical_spec, spectrum_light);

              // Simplified pre-Score, as Alpha and Beta ions are not yet tracked
              float pre_score = preScore(matched_fragments_theor_spec.size(), theoretical_spec.size());
              cout << "Number of matched peaks to theor. spectrum: " << matched_fragments_theor_spec.size() << endl;
              cout << "Number of theoretical ions: " << theoretical_spec.size() << endl;
              cout << "Pre Score: " << pre_score << endl;
              //cout << "Peptide size: " << a->second.size() << "\t" << b->second.size() << "\t" << "K Pos:" << K_pos_a[i] << "\t" << K_pos_b[i] << endl;
              if (pre_score > pScoreMax) pScoreMax = pre_score;
            }
          }
          else // enumeration mode
          {
            // determine MS2 precursors that match to the current peptide mass
            multimap<double, pair<AASequence, AASequence> >::const_iterator low_it;
            multimap<double, pair<AASequence, AASequence> >::const_iterator up_it;

            if (precursor_mass_tolerance_unit_ppm) // ppm
            {
              low_it = enumerated_cross_link_masses.lower_bound(precursor_mass - precursor_mass * precursor_mass_tolerance * 1e-6);
              up_it = enumerated_cross_link_masses.upper_bound(precursor_mass + precursor_mass * precursor_mass_tolerance * 1e-6);
            }
            else // Dalton
            {
              low_it = enumerated_cross_link_masses.lower_bound(precursor_mass - precursor_mass_tolerance);
              up_it =  enumerated_cross_link_masses.upper_bound(precursor_mass + precursor_mass_tolerance);
            }

            if (low_it == up_it) continue; // no matching precursor in data

            // loop over cross-link candidates
            for (; low_it != up_it; ++low_it)
            {
              const pair<AASequence, AASequence>& candidate = low_it->second;

//            cout << "Pair: " << a->second << ", " << b->second << " matched to light spectrum " << scan_index_light << " with m/z: " << spectrum_light.getPrecursors()[0].getMZ() <<  endl;
//            cout << a->second.getMonoWeight() << ", " << b->second.getMonoWeight() << " cross_link_mass: " <<  cross_link_mass <<  endl;


	      TheoreticalSpectrumGeneratorXLinks specGen;
	      RichPeakSpectrum theoretical_spec_beta;
	      RichPeakSpectrum theoretical_spec_alpha;
	      RichPeakSpectrum theoretical_spec_xlinks_alpha;
	      RichPeakSpectrum theoretical_spec_xlinks_beta;

              TheoreticalSpectrumGeneratorXLinks::ProteinProteinCrossLink cross_link_candidate;
              cross_link_candidate.alpha = candidate.first;
              cross_link_candidate.beta = candidate.second;
              cross_link_candidate.cross_link_position.first = candidate.first.toString().find('K');
              cross_link_candidate.cross_link_position.second = candidate.second.toString().find('K');

              //specGen.getSpectrum(theoretical_spec, cross_link_candidate, 1);
              //getXLinkIonSpectrum(theoretical_spec_xlinks , cross_link_candidate, 1)
              specGen.getCommonIonSpectrum(theoretical_spec_alpha, cross_link_candidate.alpha, 3);
              specGen.getCommonIonSpectrum(theoretical_spec_beta, cross_link_candidate.beta, 3);
              specGen.getXLinkIonSpectrum(theoretical_spec_xlinks_alpha, theoretical_spec_xlinks_beta, cross_link_candidate, 4, 6);

              std::vector< std::pair< Size, Size > > matched_spec_alpha;
              std::vector< std::pair< Size, Size > > matched_spec_beta;
              std::vector< std::pair< Size, Size > > matched_spec_xlinks_alpha;
              std::vector< std::pair< Size, Size > > matched_spec_xlinks_beta;

              ms2_aligner.getSpectrumAlignment(matched_spec_alpha, theoretical_spec_alpha, ps.spectra_common_peaks[scan_index]);
              ms2_aligner.getSpectrumAlignment(matched_spec_beta, theoretical_spec_beta, ps.spectra_common_peaks[scan_index]);
              ms2_aligner.getSpectrumAlignment(matched_spec_xlinks_alpha, theoretical_spec_xlinks_alpha, ps.spectra_xlink_peaks[scan_index]);
              ms2_aligner.getSpectrumAlignment(matched_spec_xlinks_beta, theoretical_spec_xlinks_beta, ps.spectra_xlink_peaks[scan_index]);

              /*
              ms2_aligner.getSpectrumAlignment(matched_spec_alpha, theoretical_spec_alpha, spectrum_light);
              ms2_aligner.getSpectrumAlignment(matched_spec_beta, theoretical_spec_beta, spectrum_light);
              ms2_aligner.getSpectrumAlignment(matched_spec_xlinks_alpha, theoretical_spec_xlinks_alpha, spectrum_light);
              ms2_aligner.getSpectrumAlignment(matched_spec_xlinks_beta, theoretical_spec_xlinks_beta, spectrum_light);
              */

              // Pre-Score calculations
              Size matched_alpha_count = matched_spec_alpha.size() + matched_spec_xlinks_alpha.size();
              Size theor_alpha_count = theoretical_spec_alpha.size() + theoretical_spec_xlinks_alpha.size();
              Size matched_beta_count = matched_spec_beta.size() + matched_spec_xlinks_beta.size();
              Size theor_beta_count = theoretical_spec_beta.size() + theoretical_spec_xlinks_beta.size();

              if (matched_alpha_count + matched_beta_count > 0)
              {
                  // Simplified pre-Score
                  //float pre_score = preScore(matched_fragments_theor_spec.size(), theoretical_spec.size());
                  float pre_score = preScore(matched_alpha_count, theor_alpha_count, matched_beta_count, theor_beta_count);
                  cout << "Number of matched peaks to theor. spectrum: " << matched_alpha_count << "\t" << matched_beta_count << endl;
                  cout << "Number of theoretical ions: " << theor_alpha_count << "\t" << theor_beta_count << endl;
                  cout << "Pre Score: " << pre_score << endl;
                  //cout << "Peptide size: " << a->second.size() << "\t" << b->second.size() << "\t" << "K Pos:" << K_pos_a[i] << "\t" << K_pos_b[i] << endl;
                  if (pre_score > pScoreMax) pScoreMax = pre_score;

                  // %TIC score calculations
                  double matched_current = 0;
                  for (SignedSize i = 0; i < (SignedSize)matched_spec_alpha.size(); ++i) matched_current += ps.spectra_common_peaks[scan_index][matched_spec_alpha[i].second].getIntensity();
                  for (SignedSize i = 0; i < (SignedSize)matched_spec_beta.size(); ++i) matched_current += ps.spectra_common_peaks[scan_index][matched_spec_beta[i].second].getIntensity();
                  for (SignedSize i = 0; i < (SignedSize)matched_spec_xlinks_alpha.size(); ++i) matched_current += ps.spectra_xlink_peaks[scan_index][matched_spec_xlinks_alpha[i].second].getIntensity();
                  for (SignedSize i = 0; i < (SignedSize)matched_spec_xlinks_beta.size(); ++i) matched_current += ps.spectra_xlink_peaks[scan_index][matched_spec_xlinks_beta[i].second].getIntensity();

                  double total_current = 0;
                  for (SignedSize i = 0; i < (SignedSize)spectrum_light.size(); ++i) total_current += spectrum_light[i].getIntensity();

                  double TIC = matched_current / total_current * 100.0;
                  cout << "matched current: " << matched_current << "\t total_current: " << total_current << endl;
                  cout << "%TIC= " << TIC << "%" << endl;
                  if (TIC > TICMax) TICMax = TIC;

                  // match-odds score
                  double range_c_alpha = theoretical_spec_alpha[theoretical_spec_alpha.size()-1].getMZ() -  theoretical_spec_alpha[0].getMZ();
                  double range_x_alpha= theoretical_spec_xlinks_alpha[theoretical_spec_xlinks_alpha.size()-1].getMZ() -  theoretical_spec_xlinks_alpha[0].getMZ();
                  double range_c_beta = theoretical_spec_beta[theoretical_spec_beta.size()-1].getMZ() -  theoretical_spec_beta[0].getMZ();
                  double range_x_beta = theoretical_spec_xlinks_beta[theoretical_spec_xlinks_beta.size()-1].getMZ() -  theoretical_spec_xlinks_beta[0].getMZ();
                  double a_priori_ca = aPrioriProb(fragment_mass_tolerance, theoretical_spec_alpha.size(), range_c_alpha, 3);
                  double a_priori_xa = aPrioriProb(fragment_mass_tolerance, theoretical_spec_xlinks_alpha.size(), range_x_alpha, 3);
                  double a_priori_cb = aPrioriProb(fragment_mass_tolerance, theoretical_spec_beta.size(), range_c_beta, 3);
                  double a_priori_xb = aPrioriProb(fragment_mass_tolerance, theoretical_spec_xlinks_beta.size(), range_x_beta, 3);
                  double match_odds_c_alpha = -log(1 - cumulativeBinomial(theoretical_spec_alpha.size(), matched_spec_alpha.size(), a_priori_ca) );
                  double match_odds_x_alpha = -log(1 - cumulativeBinomial(theoretical_spec_xlinks_alpha.size(), matched_spec_xlinks_alpha.size(), a_priori_xa) );
                  double match_odds_c_beta = -log(1 - cumulativeBinomial(theoretical_spec_beta.size(), matched_spec_beta.size(), a_priori_cb) );
                  double match_odds_x_beta = -log(1 - cumulativeBinomial(theoretical_spec_xlinks_beta.size(), matched_spec_xlinks_beta.size(), a_priori_xb) );

                  double match_odds = (match_odds_c_alpha + match_odds_x_alpha + match_odds_c_beta + match_odds_x_beta) / 4;

                  cout << "Range Alpha: " << theoretical_spec_alpha[theoretical_spec_alpha.size()-1].getMZ()  << " - " << (double) theoretical_spec_alpha[0].getMZ() << " = " << range_c_alpha << "\t A Priori probaility common alpha: " << a_priori_ca << endl;
                  cout << "Match-Odds Score Alpha: " << match_odds_c_alpha << "\t Match-Odds Final: " << match_odds << endl;
                  if (match_odds > matchOddsMax) matchOddsMax = match_odds;
              }
            }              
          }


          // TODO generate theoretical spectrum for cross link ions (if possible and neccessary)
          // TODO pre scoring done only with common peaks
          // TODO Alpha / Beta designation based on order, should be based on length / mass of peptide ?
          // TODO Charge state is a constant = 1, highest pre_scores without xlink peaks
          // TODO speed up, alpha chain spectrum can be computed one time for all betas...
          // TODO for each precursor peak pair in MS1, find 50 a+b combos with highest pre_score

          /**  Old Spectrumgenerator
          TheoreticalSpectrumGenerator specGen;

          // Common Peaks with old Spectrum Generator
          RichPeakSpectrum theoretical_alpha_chain;
          RichPeakSpectrum theoretical_beta_chain;
          specGen.getSpectrum(theoretical_alpha_chain, a->second, 1);
          specGen.getSpectrum(theoretical_beta_chain, b->second, 1);

          std::vector< std::pair< Size, Size > > matched_fragments_alpha_chain;
          std::vector< std::pair< Size, Size > > matched_fragments_beta_chain;
          ms2_aligner.getSpectrumAlignment(matched_fragments_alpha_chain, theoretical_alpha_chain, spectrum_light);
          ms2_aligner.getSpectrumAlignment(matched_fragments_beta_chain, theoretical_beta_chain, spectrum_light);

          float pre_score = preScore(matched_fragments_alpha_chain.size(), theoretical_alpha_chain.size(), matched_fragments_beta_chain.size(), theoretical_beta_chain.size());
          cout << "Numbers of matched peaks to theor. spectra: " << matched_fragments_alpha_chain.size() << "\t" << matched_fragments_beta_chain.size() << endl;
          cout << "Numbers of theoretical ions: " << theoretical_alpha_chain.size() << "\t" << theoretical_beta_chain.size() << endl;
          cout << "Pre Score: " << pre_score << endl;
          cout << "Peptide size: " << a->second.size() << "\t" << b->second.size() << "\t" << "K Pos:" << a->second.toString().find("K") << "\t" << b->second.toString().find("K") << endl;
          if (pre_score > pScoreMax) pScoreMax = pre_score;
          **/
	  /**
          for (Size i = 0; i < K_pos_a.size(); ++i)
          {
            for (Size j = 0; j < K_pos_b.size(); ++j)
            {

              // New TheoreticalSpectrumGeneratorXLinks
              // Common and XLink Peaks with new Spectrum Generator TheoreticalSpectrumGeneratorXLinks
              TheoreticalSpectrumGeneratorXLinks specGen;
              RichPeakSpectrum theoretical_spec;
              specGen.getSpectrum(theoretical_spec, a->second, b->second, K_pos_a[i], K_pos_b[i], 1);

              std::vector< std::pair< Size, Size > > matched_fragments_theor_spec;
              ms2_aligner.getSpectrumAlignment(matched_fragments_theor_spec, theoretical_spec, spectrum_light);

              // Simplified pre-Score, as Alpha and Beta ions are not yet tracked
              float pre_score = preScore(matched_fragments_theor_spec.size(), theoretical_spec.size());
              cout << "Number of matched peaks to theor. spectrum: " << matched_fragments_theor_spec.size() << endl;
              cout << "Number of theoretical ions: " << theoretical_spec.size() << endl;
              cout << "Pre Score: " << pre_score << endl;
              cout << "Peptide size: " << a->second.size() << "\t" << b->second.size() << "\t" << "K Pos:" << K_pos_a[i] << "\t" << K_pos_b[i] << endl;
              if (pre_score > pScoreMax) pScoreMax = pre_score;
            }
          }
	  **/
        }
      }

    }
/**

    // calculate mass pairs
    for (map<StringView, AASequence>::const_iterator a = processed_peptides.begin(); a != processed_peptides.end(); ++a)
    {
      if (++counter % 1000 == 0) cout << counter * 100.0 / processed_peptides.size() << endl;

      if (a->second.toString().find("K") >= a->second.size()-1)
      {
        continue;
      }

      for (map<StringView, AASequence>::const_iterator b = a; b != processed_peptides.end(); ++b)
      {

        if (b->second.toString().find("K") >= b->second.size()-1)
        {
          continue;
        }

        // Find all positions of lysine (K) in the peptides (possible scross-linking sites)
        vector <Size> K_pos_a;
        vector <Size> K_pos_b;
        for (Size i = 0; i < a->second.size()-1; ++i)
        {
          if (a->second.toString()[i] == 'K') K_pos_a.push_back(i);
        }
        for (Size i = 0; i < b->second.size()-1; ++i)
        {
          if (b->second.toString()[i] == 'K') K_pos_b.push_back(i);
        }

        // mass peptide1 + mass peptide2 + cross linker mass - cross link loss
        double cross_link_mass = a->second.getMonoWeight() + b->second.getMonoWeight() + cross_link_mass_light - cross_link_mass_loss_type2;

        // determine MS2 precursors that match to the current peptide mass
        multimap<double, Size>::const_iterator low_it;
        multimap<double, Size>::const_iterator up_it;

        if (precursor_mass_tolerance_unit_ppm) // ppm
        {
          low_it = multimap_mass_2_scan_index.lower_bound(cross_link_mass - cross_link_mass * precursor_mass_tolerance * 1e-6);
          up_it = multimap_mass_2_scan_index.upper_bound(cross_link_mass + cross_link_mass * precursor_mass_tolerance * 1e-6);
        }
        else // Dalton
        {
          low_it = multimap_mass_2_scan_index.lower_bound(cross_link_mass - precursor_mass_tolerance);
          up_it = multimap_mass_2_scan_index.upper_bound(cross_link_mass + precursor_mass_tolerance);
        }

        if (low_it == up_it) continue; // no matching precursor in data

        for (; low_it != up_it; ++low_it)
        {
          const Size scan_index_light = low_it->second;
          const PeakSpectrum& spectrum_light = spectra[scan_index_light];

//          cout << "Pair: " << a->second << ", " << b->second << " matched to light spectrum " << scan_index_light << " with m/z: " << spectrum_light.getPrecursors()[0].getMZ() <<  endl;
//          cout << a->second.getMonoWeight() << ", " << b->second.getMonoWeight() << " cross_link_mass: " <<  cross_link_mass <<  endl;

        }  
      }     
    }
**/
    cout << "Pre Score maximum: " << pScoreMax << "\t TIC maximum: " << TICMax << "\t Match-Odds maximum: " << matchOddsMax << endl;
 
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPxQuest tool;
  return tool.main(argc, argv);
}

